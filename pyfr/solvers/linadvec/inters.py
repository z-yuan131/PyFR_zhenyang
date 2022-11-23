# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import (BaseAdvectionIntInters,
                                    BaseAdvectionMPIInters,
                                    BaseAdvectionBCInters)

class TplargsMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        rsolver = self.cfg.get('solver-interfaces', 'riemann-solver')
        #visc_corr = self.cfg.get('solver', 'viscosity-correction')
        #shock_capturing = self.cfg.get('solver', 'shock-capturing')
        self._tplargs = dict(ndims=self.ndims, nvars=self.nvars,
                             rsolver=rsolver,  c=self.c)


class LinearAdvectionIntInters(TplargsMixin, BaseAdvectionIntInters):
    def __init__(self, be, lhs, rhs, elemap, cfg):
        super().__init__(be, lhs, rhs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.linadvec.kernels.intconu')

        # Generate the left and right hand side view matrices
        self._base_scal_lhs = self._scal_view(lhs, 'get_base_scal_fpts_for_inter')
        self._base_scal_rhs = self._scal_view(rhs, 'get_base_scal_fpts_for_inter')

        # Generate the additional view matrices
        self._base_vect_lhs = self._vect_view(lhs, 'get_base_vect_fpts_for_inter')
        self._base_vect_rhs = self._vect_view(rhs, 'get_base_vect_fpts_for_inter')

        # Additional kernel constants
        self.c |= cfg.items_as('solver-interfaces', float)

        # Common solution at our inner interface
        self.kernels['con_u'] = lambda: be.kernel(
            'intconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._base_scal_lhs, urin=self._base_scal_rhs,
            ulout=self._base_vect_lhs, urout=self._base_vect_rhs
        )

    """
    def _gen_perm(self, lhs, rhs):
        # In the special case of β = -0.5 it is better to sort by the
        # RHS interface; otherwise we simply opt for the LHS
        beta = self.cfg.getfloat('solver-interfaces', 'ldg-beta')
        side = lhs if beta != -0.5 else rhs

        # Compute the relevant permutation
        self._perm = self._get_perm_for_view(side, 'get_base_scal_fpts_for_inter')
    """

class LinearAdvectionMPIInters(TplargsMixin, BaseAdvectionMPIInters):
    def __init__(self, be, lhs, rhsrank, rallocs, elemap, cfg):
        super().__init__(be, lhs, rhsrank, rallocs, elemap, cfg)

        be.pointwise.register('pyfr.solvers.linadvec.kernels.mpiconu')

        lhsprank = rallocs.prank
        rhsprank = rallocs.mprankmap[rhsrank]

        # Generate second set of view matrices
        self._base_vect_lhs = self._vect_xchg_view(lhs, 'get_base_vect_fpts_for_inter')
        self._base_vect_rhs = be.xchg_matrix_for_view(self._vect_lhs)

        # Additional kernel constants
        self.c |= cfg.items_as('solver-interfaces', float)

        # We require cflux(l,r,n_l) = -cflux(r,l,n_r) and
        # conu(l,r) = conu(r,l) and where l and r are left and right
        # solutions at an interface and n_[l,r] are physical normals.
        # The simplest way to enforce this at an MPI interface is for
        # one side to take β = -β for the cflux and conu kernels. We
        # pick this side (arbitrarily) by comparing the physical ranks
        # of the two partitions.
        if (lhsprank + rhsprank) % 2:
            self.c['ldg-beta'] *= 1.0 if lhsprank > rhsprank else -1.0
        else:
            self.c['ldg-beta'] *= 1.0 if rhsprank > lhsprank else -1.0

        # Allocate a tag
        vect_fpts_tag = next(self._mpi_tag_counter)

        # If we need to send our gradients to the RHS
        if self.c['ldg-beta'] != -0.5:
            self.kernels['vect_fpts_pack'] = lambda: be.kernel(
                'pack', self._base_vect_lhs
            )
            self.mpireqs['vect_fpts_send'] = lambda: self._base_vect_lhs.sendreq(
                self._rhsrank, vect_fpts_tag
            )

        # If we need to recv gradients from the RHS
        if self.c['ldg-beta'] != 0.5:
            self.mpireqs['vect_fpts_recv'] = lambda: self._base_vect_rhs.recvreq(
                self._rhsrank, vect_fpts_tag
            )
            self.kernels['vect_fpts_unpack'] = lambda: be.kernel(
                'unpack', self._base_vect_rhs
            )

        # Common solution at our MPI interface
        self.kernels['con_u'] = lambda: be.kernel(
            'mpiconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            ulin=self._base_scal_lhs, urin=self._base_scal_rhs,
            ulout=self._base_vect_lhs
        )


class LinearAdvectionBCInters(TplargsMixin, BaseAdvectionBCInters):
    def __init__(self, be, lhs, elemap, cfgsect, cfg):
        super().__init__(be, lhs, elemap, cfgsect, cfg)

        be.pointwise.register('pyfr.solvers.navstokes.kernels.bcconu')

        # Additional BC specific template arguments
        self._tplargs['bctype'] = self.type

        # Generate the left and right hand side view matrices
        self._base_scal_lhs = self._scal_view(lhs, 'get_base_scal_fpts_for_inter')

        # Additional view matrices
        self._base_vect_lhs = self._vect_view(lhs, 'get_base_vect_fpts_for_inter')

        # Additional kernel constants
        self.c |= cfg.items_as('solver-interfaces', float)

        # Common solution at our boundary interface
        self.kernels['con_u'] = lambda: be.kernel(
            'bcconu', tplargs=self._tplargs, dims=[self.ninterfpts],
            extrns=self._external_args, ulin=self._base_scal_lhs,
            ulout=self._base_vect_lhs, nlin=self._norm_pnorm_lhs,
            **self._external_vals
        )
