import numpy as np

from pyfr.solvers.baseadvecdiff import (BaseAdvectionDiffusionBCInters,
                                        BaseAdvectionDiffusionIntInters,
                                        BaseAdvectionDiffusionMPIInters)
from pyfr.solvers.euler.inters import (FluidIntIntersMixin,
                                       FluidMPIIntersMixin,
                                       MassFlowBCMixin)
from pyfr.solvers.navstokes.inters import TplargsMixin


class BaseflowIntInters:
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        # Baseflow state and gradients at interface flux points.
        self._baseflow_lhs = self._scal_view(lhs, 'get_baseflow_for_inter')
        self._baseflow_rhs = self._scal_view(rhs, 'get_baseflow_for_inter')
        self._basegrad_lhs = self._vect_view(lhs, 'get_basegrad_for_inter')
        self._basegrad_rhs = self._vect_view(rhs, 'get_basegrad_for_inter')

        """
        self._set_external('ubl', f'in fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('ubr', f'in fpdtype_t[{self.nvars}]',
                           value=self._baseflow_rhs)
        self._set_external('gbl', f'in fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)
        self._set_external('gbr', f'in fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_rhs)
        """
        self._set_external('ub', f'in view fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('gb', f'in view fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)


class BaseflowMPIInters:
    def __init__(self, be, lhs, rhsrank, elemap, cfg):
        super().__init__(be, lhs, rhsrank, elemap, cfg)

        # Prototype assumption: baseflow is smooth; use lhs for rhs.
        self._baseflow_lhs = self._scal_xchg_view(lhs, 'get_baseflow_for_inter')
        self._baseflow_rhs = self._baseflow_lhs
        self._basegrad_lhs = self._vect_xchg_view(lhs, 'get_basegrad_for_inter')
        self._basegrad_rhs = self._basegrad_lhs

        """
        self._set_external('ubl', f'in fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('ubr', f'in fpdtype_t[{self.nvars}]',
                           value=self._baseflow_rhs)
        self._set_external('gbl', f'in fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)
        self._set_external('gbr', f'in fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_rhs)
        """
        self._set_external('ub', f'in view fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('gb', f'in view fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)

        self._be.pointwise.register('pyfr.solvers.lnavstokes.kernels.bcconub')

        self.kernels['con_ub'] = lambda: self._be.kernel(
            'bcconub', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )


class BaseflowBCInters:
    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self._baseflow_lhs = self._scal_view(lhs, 'get_baseflow_for_inter')
        self._basegrad_lhs = self._vect_view(lhs, 'get_basegrad_for_inter')

        """
        self._set_external('ubl', f'in fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('gbl', f'in fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)
        """
        self._set_external('ub', f'in view fpdtype_t[{self.nvars}]',
                           value=self._baseflow_lhs)
        self._set_external('gb', f'in view fpdtype_t[{self.ndims}][{self.nvars}]',
                           value=self._basegrad_lhs)



class LinearNavierStokesIntInters(TplargsMixin,
                                  BaseflowIntInters,
                                  FluidIntIntersMixin,
                                  BaseAdvectionDiffusionIntInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.intconu')
        self._be.pointwise.register('pyfr.solvers.lnavstokes.kernels.intcflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs,
            ulout=self._comm_lhs, urout=self._comm_rhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'intcflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args,
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            nl=self._pnorm_lhs, **self._external_vals
        )


class LinearNavierStokesMPIInters(TplargsMixin,
                                  BaseflowMPIInters,
                                  FluidMPIIntersMixin,
                                  BaseAdvectionDiffusionMPIInters):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._be.pointwise.register('pyfr.solvers.navstokes.kernels.mpiconu')
        self._be.pointwise.register('pyfr.solvers.lnavstokes.kernels.mpicflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._scal_lhs, urin=self._scal_rhs, ulout=self._comm_lhs
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'mpicflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args,
            ul=self._scal_lhs, ur=self._scal_rhs,
            gradul=self._vect_lhs, gradur=self._vect_rhs,
            artviscl=self._artvisc_lhs, artviscr=self._artvisc_rhs,
            nl=self._pnorm_lhs, **self._external_vals
        )


class LinearNavierStokesBaseBCInters(TplargsMixin, 
                                     BaseflowBCInters,
                                     BaseAdvectionDiffusionBCInters):
    cflux_state = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type
        self._tplargs['bccfluxstate'] = self.cflux_state

        self._be.pointwise.register('pyfr.solvers.lnavstokes.kernels.bcconu')
        self._be.pointwise.register('pyfr.solvers.lnavstokes.kernels.bccflux')

        self.kernels['con_u'] = lambda: self._be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._scal_lhs,
            ulout=self._comm_lhs, nlin=self._pnorm_lhs,
            **self._external_vals
        )
        self.kernels['comm_flux'] = lambda: self._be.kernel(
            'bccflux', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ul=self._scal_lhs,
            gradul=self._vect_lhs, nl=self._pnorm_lhs,
            artviscl=self._artvisc_lhs, **self._external_vals
        )

        if self._ef_enabled:
            self._be.pointwise.register(
                'pyfr.solvers.navstokes.kernels.bccent'
            )

            self.kernels['comm_entropy'] = lambda: self._be.kernel(
                'bccent', tplargs=self._tplargs, dims=[self.ninterfpts],
                extrns=self._external_args, entmin_lhs=self._entmin_lhs,
                nl=self._pnorm_lhs, ul=self._scal_lhs, **self._external_vals
            )


class LinearNavierStokesNoSlpIsotWallBCInters(LinearNavierStokesBaseBCInters):
    type = 'no-slp-isot-wall'
    cflux_state = 'ghost-imperm'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c['cpTw'], = self._eval_opts(['cpTw'])
        self.c |= self._exp_opts('uvw'[:self.ndims], lhs,
                                 default={'u': 0, 'v': 0, 'w': 0})


class LinearNavierStokesNoSlpAdiaWallBCInters(LinearNavierStokesBaseBCInters):
    type = 'no-slp-adia-wall'
    cflux_state = 'ghost-imperm'


class LinearNavierStokesSlpAdiaWallBCInters(LinearNavierStokesBaseBCInters):
    type = 'slp-adia-wall'
    cflux_state = None


class LinearNavierStokesCharRiemInvBCInters(LinearNavierStokesBaseBCInters):
    type = 'char-riem-inv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class LinearNavierStokesSupInflowBCInters(LinearNavierStokesBaseBCInters):
    type = 'sup-in-fa'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'p', 'u', 'v', 'w'][:self.ndims + 2], lhs
        )


class LinearNavierStokesSupOutflowBCInters(LinearNavierStokesBaseBCInters):
    type = 'sup-out-fn'
    cflux_state = 'ghost'


class LinearNavierStokesSubInflowFrvBCInters(LinearNavierStokesBaseBCInters):
    type = 'sub-in-frv'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(
            ['rho', 'u', 'v', 'w'][:self.ndims + 1], lhs,
            default={'u': 0, 'v': 0, 'w': 0}
        )


class LinearNavierStokesSubInflowFtpttangBCInters(LinearNavierStokesBaseBCInters):
    type = 'sub-in-ftpttang'
    cflux_state = 'ghost'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        gamma = self.cfg.getfloat('constants', 'gamma')

        # Pass boundary constants to the backend
        self.c['cpTt'], = self._eval_opts(['cpTt'])
        self.c['pt'], = self._eval_opts(['pt'])
        self.c['Rdcp'] = (gamma - 1.0)/gamma

        # Calculate u, v velocity components from the inflow angle
        theta = self._eval_opts(['theta'])[0]*np.pi/180.0
        velcomps = np.array([np.cos(theta), np.sin(theta), 1.0])

        # Adjust u, v and calculate w velocity components for 3-D
        if self.ndims == 3:
            phi = self._eval_opts(['phi'])[0]*np.pi/180.0
            velcomps[:2] *= np.sin(phi)
            velcomps[2] *= np.cos(phi)

        self.c['vc'] = velcomps[:self.ndims]


class LinearNavierStokesSubOutflowBCInters(LinearNavierStokesBaseBCInters):
    type = 'sub-out-fp'
    cflux_state = 'ghost'

    def __init__(self, be, lhs, elemap, cfgsect, cfg, bccomm):
        super().__init__(be, lhs, elemap, cfgsect, cfg, bccomm)

        self.c |= self._exp_opts(['p'], lhs)


class LinearNavierStokesCharRiemInvMassFlowBCInters(MassFlowBCMixin,
                                              LinearNavierStokesBaseBCInters):
    type = 'char-riem-inv-mass-flow'
    cflux_state = 'ghost'
