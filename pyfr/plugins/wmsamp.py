from collections import namedtuple
import numpy as np
import os

from pyfr.mpiutil import get_comm_rank_root, mpi
from pyfr.plugins.base import BaseSolverPlugin, init_csv
from pyfr.points import PointLocator, PointSampler
from pyfr.quadrules import get_quadrule
from pyfr.plugins.sampler import _process_con_to_pri
from pyfr.nputil import npeval



class WMLESSamplingHeightController:
    def __init__(self, ym_init, #y1,
                 yplus_target=50.0,
                 Tf=10.0,        # filter time scale
                 beta=0.02):     # adaptation relaxation
        #self.ym = ym_init
        #self.y1 = y1

        self.yplus_target = yplus_target
        self.Tf = Tf
        self.beta = beta

        self.utau_filt = None

        # Hard bounds
        #self.ym_min = 1.0 * y1
        #self.ym_max = 3.0 * y1

    def update(self, utau_inst, nu, dt, ym):
        eps = 1e-12

        # --- Low-pass filter u_tau ---
        if self.utau_filt is None:
            self.utau_filt = utau_inst
        else:
            alpha = dt / (self.Tf + dt)
            self.utau_filt = (
                (1.0 - alpha) * self.utau_filt
                + alpha * utau_inst
            )

        # --- Desired ym from filtered u_tau ---
        ym_star = self.yplus_target * nu / (self.utau_filt + eps)

        # --- Relaxation update ---
        ym = (
            (1.0 - self.beta) * ym
            + self.beta * ym_star
        )

        # --- Hard bounds ---
        #self.ym = max(self.ym_min, min(self.ym, self.ym_max))

        return ym

    def current_yplus(self, ym, nu):
        return ym * self.utau_filt / (nu + 1e-12)



class WMSampPlugin(BaseSolverPlugin):
    name = 'wmsamp'
    systems = ['*']
    formulations = ['std']
    dimensions = [2, 3]

    def __init__(self, intg, cfgsect):
        super().__init__(intg, cfgsect)

        # Underlying elements class
        self.elementscls = intg.system.elementscls

        # Also element maps
        self.ele_map = intg.system.ele_map

        # Data type
        self.dtype = intg.system.backend.fpdtype

        # Get sampler for the edge locations
        self._get_wmloc(intg, cfgsect)

        # If we get the stress 
        self.adym = self.cfg.getbool(cfgsect, 'adapt-sample', False)
        if self.adym:

            self.utau_avg = 0
            
            ym = self.cfg.getfloat('constants', 'ym')
            Tf = self.cfg.getfloat(cfgsect, 'filter-timescale')
            self.ymctl = WMLESSamplingHeightController(ym , Tf)

            self.dt_avg = self.cfg.getfloat(cfgsect, 'dt-avg')
            #self.tupdate_next = intg.tcurr + self.dt_avg
            self.tupdate_next = intg.tcurr 

        #"""
        # If we get the stress 
        self.getstress = self.cfg.getbool(cfgsect, 'get-stress', False)
        if self.getstress:
            self.dt_out = self.cfg.getfloat(cfgsect, 'dt-out')
            self.dt_samp_out = self.cfg.getfloat(cfgsect, 'dt-sample-out')
            self.tout_last = intg.tcurr

            self.dt_samp_out += intg.tcurr
            self.dt_freq = self.dt_samp_out / self.dt_out

            if not intg.isrestart:
                self.tout_last -= self.dt_out


            # Root rank should open the file
            comm, rank, root = get_comm_rank_root()
            if rank == root:
                self._init_csv(intg)
                # write the coordinate
                nploc = np.concatenate(self.nploc, axis = 0)
                for k in range(self.ndims):
                    self.csv(intg.tcurr, *nploc[:,k])

                # avg of the stress mag 
                self.tt = 0
        #"""
            
        
        # Prepare WM kernels
        self.WMInfo, self.stressInfo = self._wmbc_kern(intg)

    def __call__(self, intg):
        # For debug
        comm, rank, root = get_comm_rank_root()

        print(rank, 'a', flush = True)


        # Updates force and locations
        nsoln = self._update_soln(intg)

        # Updates
        for wminfo, soln in zip(self.WMInfo, nsoln):
            wminfo.set(soln.T)

        """
        # If adpative stress
        if self.adym:
            for stinfo in self.stressInfo:
                info = stinfo.get()
                ym, utau, nu = info[0], info[1], info[2]

                #print(np.max(info), np.min(info))
                #return
                #raise RuntimeError

                if intg.tcurr - self.tupdate_next > self.tol:
                    dt = intg.tcurr - self.tupdate_next
                    ymn = self.ymctl.update(utau, nu, dt, ym)

                    print(np.min(ym),np.max(self.ymctl.current_yplus(ymn, nu)),np.min(self.ymctl.current_yplus(ymn, nu)))

                    # Updates
                    info[0] = ymn 
                    stinfo.set(info)

                    self.tupdate_next = intg.tcurr 
        """


        """
                if intg.tcurr < self.tupdate_next:
                    # Accumulate 
                    self.utau_avg += utau/nu 

                else:
                    # Update ym
                    ymn = self._update_ym(self.utau_avg/self.dt_avg, ym)
                    info[0] = ymn 
                    stinfo.set(info)

                    self.tupdate_next = intg.tcurr + self.dt_avg
                    self.utau_avg = 0
        """

        # Output the stress information 
        if self.getstress:
            if intg.tcurr - self.tout_last < self.dt_out - self.tol:
                return 

            utinfo = []
            for _id, stinfo in enumerate(self.stressInfo):
                # Gather the info from the ranks with the wall model bc
                bccomm = self.wmbc_comms[_id]
                info = stinfo.get()
                info = bccomm.gather(info, root=root)

                if rank == root:
                    info = np.concatenate(info, axis = 1)

                    #for k in tauw:
                    #    self.csv(intg.tcurr, *k)
                    self.tt += info

                    #print('avg', intg.tcurr, self.dt_samp_out)

                    if intg.tcurr > self.dt_samp_out - self.tol:
                        print('rank 0 flushing')
                        self.tt = self.tt / self.dt_freq
                        for k in [1, 2]:
                            self.csv(intg.tcurr, *self.tt[k])

                        self.dt_samp_out += intg.tcurr
                        self.tt = 0


            # Update the last output time
            self.tout_last = intg.tcurr



    def _init_csv(self, intg):
        self.csv = init_csv(self.cfg, self.cfgsect, self._header(intg),
                            nflush=len(self.nploc))

    def _header(self, intg):
        dims = 'xyz'[:self.ndims]

        vmap = ['tau0', 'tau1', 'tau2'][:self.ndims]

        colnames = ['t', *dims, *vmap]

        return ','.join(colnames)


    def _update_soln(self, intg):
        # Fetch the solution
        soln = list(intg.soln)

        # MPI
        comm, rank, root = get_comm_rank_root()

        sln = []
        for psampler, bcidx in zip(self.psampler, self.bcidx):
            nsoln = psampler.sample(soln, process=self._process)

            # Broadcast interolated solution to all ranks
            nsoln = comm.bcast(nsoln, root = root)

            # Locate to the boundary points within the rank
            start, end = bcidx[rank]
            sln.append(nsoln[start:end])
        return sln

    def _wmbc_kern(self, intg):
        wmbcm, tauw_exp = [], []
        for npts, inters, ym in zip(self.npts, self.inters, self.ym): 
            if not inters:
                continue

            # Initialize data to broadcast
            npts, nvars = npts, self.nvars
            nsoln = np.zeros((nvars, npts))

            # Add a soln matrix to backend
            tmp = intg.backend.matrix(nsoln.shape, nsoln, tags={'align'})
            wmbcm.append(tmp)
            inters._set_external('nsoln', f'in fpdtype_t[{nvars}]', value=tmp)

            # Add a location matrix ym to backend
            _ym = ym.reshape(1, -1)
            #a = np.zeros()
            _ym = intg.backend.matrix(_ym.shape, _ym, tags={'align'})
            inters._set_external('ym', f'in fpdtype_t', value=_ym)

            # If we want to output the stress 
            if self.adym or self.getstress:
                tau = np.zeros((3, npts))

                tauw = intg.backend.matrix(tau.shape, tau, tags={'align'})
                tauw_exp.append(tauw)
                inters._set_external('tauw_exp', f'inout fpdtype_t[3]', value=tauw)

        return wmbcm, tauw_exp

    def _get_wmloc(self, intg, cfgsect):
        # Load the boundary layer element 
        bcon = intg.system.mesh.bcon 

        # See which ranks have the boundary
        comm, rank, root = get_comm_rank_root()
        wmbcs = self.cfg.get(cfgsect, 'bc-names')
        wmbcs = [s.strip() for s in wmbcs.split(',')]
        bcranks = [suffix in bcon for suffix in wmbcs]

        # Split the communicators for ranks with the boundary
        self.wmbc_comms = wmbc_comms = []
        for _id, bcrank in enumerate(bcranks):
            color = _id if bcrank else mpi.UNDEFINED
            wmbc_comms.append(comm.Split(color=color, key=rank))

        # Interfaces for the wall model
        _bc_inters = intg.system._bc_inters
        
        self.inters, pts, self.ym = [], [], []
        its_names = {its.name:its for its in _bc_inters}
        for wmbc in wmbcs:
            if wmbc in its_names:
                # Get interface
                inters = its_names[wmbc]
                self.inters.append(inters)

                # Load ploc and normals
                lhs = intg.system.mesh.bcon[inters.name]
                ploc = inters._const_mat(lhs, 'get_ploc_for_inter').get()
                nl = inters._const_mat(lhs, 'get_pnorms_for_inter').get()

                # Normalise nl and let nl points into the element
                mag_nl = np.linalg.norm(nl, axis = 0)
                nl = -nl / mag_nl[None, :]

                # Get the location which is ym away from the wall
                ym = self._get_ym(cfgsect, ploc)
                self.ym.append(ym)
                nloc = ploc + ym * nl 
                nloc = nloc.T
            else:
                self.inters.append(None)
                self.ym.append(None)
                nloc = np.empty((0, self.ndims), dtype=self.dtype)

            pts.append(nloc)

        # Prepare global interpolation point sets for each requested BC
        self.nploc, self.bcidx, self.npts = [], [], []
        for nloc in pts:
            gnloc, gbcidx, lnpts = self._get_wmeval(nloc)
            self.nploc.append(gnloc)
            self.bcidx.append(gbcidx)
            self.npts.append(lnpts)

        """
        if rank == 0:
            import matplotlib.pyplot as plt 
            plt.figure()
            nnloc = np.concatenate(self.nploc, axis = 0)
            print(nnloc.shape)
            plt.plot(nnloc[:,0], nnloc[:,1],'.')
            plt.plot(ploc[0], ploc[1],'r.')
            plt.axis('equal')
            plt.savefig('/scratch/zhenyang/compute/pyfr/Naca0012trip/wmles/2d/wm_loc.png')
            plt.show()
        comm.barrier()
        raise RuntimeError
        #"""


        # Construct operators for the interpolation
        self.mesh = intg.system.mesh
        self._process = _process_con_to_pri(self.elementscls, self.ndims,
                                                self.cfg)
        self.psampler = []
        for nloc in self.nploc:
            psampler = PointSampler(self.mesh, nloc)
            psampler.configure_with_intg_nvars(intg, self.nvars)
            self.psampler.append(psampler)

    def _get_ym(self, cfgsect, plocs):
        # Get the sampling location ym
        vars = dict(zip('xyz', plocs))
        ym = npeval(self.cfg.getexpr(cfgsect, 'ym'), vars)
        if type(ym) is float:
            ym = ym * np.ones(plocs.shape[1], dtype=self.dtype)
        return ym

    def _get_wmeval(self, nloc):
        # Boardcast the location 
        comm, rank, root = get_comm_rank_root()

        # Local 2D array
        n_loc = np.ascontiguousarray(nloc, dtype=self.dtype)
        n_i, m = n_loc.shape[0], self.ndims

        # npts per rank
        row_counts = np.array(comm.allgather(n_i), dtype=np.int32)
        row_starts = np.zeros_like(row_counts)
        row_starts[1:] = np.cumsum(row_counts[:-1])
        row_ends = row_starts + row_counts
        
        # look-up table
        bcidx = np.column_stack((row_starts, row_ends))

        # Prepare for boardcast
        elem_counts = row_counts * m
        elem_displs = row_starts * m

        n_rows_global = int(np.sum(row_counts))
        a_global_flat = np.empty(n_rows_global * m, dtype=self.dtype)

        comm.Allgatherv(
            n_loc.ravel(),
            [a_global_flat, elem_counts, elem_displs, mpi.DOUBLE]
        )

        # Reshape to 2D global array
        nploc = a_global_flat.reshape(n_rows_global, m)
        return nploc, bcidx, n_i

    def _config_ubases(self, intg):
        # Get the solution bases from the system
        ubases = {etype: eles.basis.ubasis
                  for etype, eles in intg.system.ele_map.items()}
        
        return ubases
