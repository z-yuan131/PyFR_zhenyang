from collections import namedtuple
import numpy as np

from pyfr.plugins.solver.almcoeff import prop13x6
from pyfr.mpiutil import get_comm_rank_root
from pyfr.plugins.solver.base import BaseSolverPlugin
from pyfr.points import PointLocator, PointSampler
from pyfr.plugins.cli.sampler import _process_con_to_pri

ALMInfo = namedtuple('alminfo', ['nloc', 'forc'])


class ALMPlugin(BaseSolverPlugin):
    name = 'alm'
    systems = '.*'
    dimensions = '3'

    def __init__(self, intg, cfgsect):
        super().__init__(intg, cfgsect)

        # Underlying elements class
        self.elementscls = intg.system.elementscls

        # Also element maps
        self.ele_map = intg.system.ele_map

        # Initialize the parameters
        self._init_param(cfgsect)
        
        # Data format for ALM calculations
        self._process = _process_con_to_pri(self.elementscls, self.ndims,
                                                self.cfg)

        # Construct and configure the point sampler and locator
        self.mesh = intg.system.mesh
        self.psampler = PointSampler(self.mesh, self.pts)
        self.ubases = self._config_ubases(intg)

        self.locf = ['cidx', 'eidx', 'tloc']
        self.plocator = PointLocator(self.mesh)

        # allocate the inital time
        self.told = intg.tcurr

        # Data type
        self.dtype = intg.system.backend.fpdtype

        # Prepare ALM kernels
        self.alm_info = self._alm_kern(intg)

    def _config_ubases(self, intg):
        # Get the solution bases from the system
        ubases = {etype: eles.basis.ubasis
                  for etype, eles in intg.system.ele_map.items()}
        
        return ubases

    def _alm_kern(self, intg):
        pts = self.pts

        # Initialize data to broadcast
        self.forc = forc = np.zeros(pts.shape, self.dtype)
        nploc = np.zeros(pts.shape, self.dtype)

        npts, ndim = len(pts), self.ndims

        # Add macro kernel to backends
        alm_info = {}
        for etype, eles in self.ele_map.items():
            eles.add_src_macro(
                'pyfr.plugins.solver.kernels.alm', 'alm',
                self.macro_params, ploc=True, soln=False
            )

            alm_info[etype] = vs = ALMInfo(
                    intg.backend.matrix(pts.shape, pts, tags={'align'}),
                    intg.backend.matrix(forc.shape, forc, tags={'align'}),
                    )

            eles._set_external('forc', f'in broadcast fpdtype_t[{npts}][{ndim}]',
                               value=vs.forc)
            eles._set_external('nloc', f'in broadcast fpdtype_t[{npts}][{ndim}]', 
                               value=vs.nloc)

        return alm_info

    def _init_param(self, cfgsect):
        # Scalar ALM parameters
        e = self.cfg.getfloat(cfgsect, 'e')
        self.omega = omega = self.cfg.getfloat(cfgsect, 'omega')
        self.M = self.cfg.getfloat(cfgsect, 'M')
        self.mu = self.cfg.getfloat(cfgsect, 'mu')

        # Radial actuator-point locations, centred in each radial segment
        r_ini = self.cfg.getfloat(cfgsect, 'r_ini')
        r_end = self.cfg.getfloat(cfgsect, 'r_end')
        n_ap = int(self.cfg.getfloat(cfgsect, 'n_actuator_points'))

        n_sample = np.linspace(r_ini, r_end, n_ap + 1)
        r_sample = (n_sample[:-1] + n_sample[1:])/2
        dr = n_sample[1] - n_sample[0]

        # Repeat the radial locations for each blade
        theta_blades = np.radians(self.cfg.getliteral(cfgsect, 'thet0'))
        n_blades = len(theta_blades)
        n_pts_per_blade = len(r_sample)
        npts = n_blades * n_pts_per_blade

        self.r = r = np.tile(r_sample, n_blades)
        self.theta0 = theta0 = np.repeat(theta_blades, n_pts_per_blade)
        self.dr = np.full(npts, dr)
        self.omegar = omega*r

        # Blade pitch and chord at each actuator point
        self.betas, self.c = prop13x6.pitch(r)

        # Initial actuator point locations in the x-normal rotor plane
        self.pts = np.zeros((npts, self.ndims))
        self.pts[:, 1] = r*np.cos(theta0)
        self.pts[:, 2] = r*np.sin(theta0)

        # Constant parameters for the macro kernel
        self.macro_params = {
            'eph': e,
            'eph2': 1/e**2,
            'ephpi3': 1/e**3/(np.pi)**(3/2),
            'r': r,
            'npts': npts
        }


    def __call__(self, intg):
        # Update force
        self._update_forc(intg)
        pts, forc = self.pts, self.forc 

        for etype, info in self.alm_info.items():
            # Updates point locations and forces
            info.nloc.set(pts)
            info.forc.set(forc)

        # Reset time
        self.told = intg.tcurr
        
    def _update_solution(self, intg): 
        # New location
        theta = self.theta0 + self.omega*intg.tcurr
        self.pts[:, 1] = self.r*np.cos(theta)
        self.pts[:, 2] = self.r*np.sin(theta)

        # Locate the new point list
        locs = self.plocator.locate(self.pts)[self.locf]

        # Fetch the solution
        soln = list(intg.soln)

        # Sample the solution
        self.psampler.locs = locs 
        self.psampler._configure_ubases_nvars(self.ubases, self.nvars)
        samps = self.psampler.sample(soln, process=self._process)

        # Broadcast to all ranks
        comm, rank, root = get_comm_rank_root()
        samps = comm.bcast(samps, root = root)

        return samps
    
    def _update_forc(self, intg): 
        # Primitive flow samples at each actuator point
        nsoln = np.asarray(self._update_solution(intg))
        rho = nsoln[:, 0]
        u = nsoln[:, 1]
        v = nsoln[:, 2]
        w = nsoln[:, 3]

        theta = self.theta0 + self.omega*intg.tcurr
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        # Velocity components in the axial and tangential blade frame
        vtheta_rel = self.omegar + v*sin_theta - w*cos_theta

        # Inflow angle, effective angle of attack, and relative speed
        phi = np.arctan2(u, vtheta_rel)
        alpha = self.betas - phi

        vrel2 = u * u + vtheta_rel * vtheta_rel
        vrel = np.sqrt(vrel2)
        reynolds = self.c * rho * vrel / self.mu

        # Aerodynamic coefficients and sectional forces
        cl_raw, cd_raw = prop13x6.clcd(np.degrees(alpha), reynolds)
        compress = np.sqrt(1 - vrel2)
        lift = 0.5 * rho * vrel2 * self.c * (cl_raw / compress) * self.dr
        drag = 0.5 * rho * vrel2 * self.c * (cd_raw / compress) * self.dr

        # Transform blade-frame forces back to solver coordinates
        ftheta = drag*np.cos(phi) + lift*np.sin(phi)
        fx = -lift*np.cos(phi) + drag*np.sin(phi)
        fy = ftheta*sin_theta
        fz = -ftheta*cos_theta

        self.forc[:, 0] = -fx
        self.forc[:, 1] = -fy
        self.forc[:, 2] = -fz
        
