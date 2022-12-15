# -*- coding: utf-8 -*-

from pyfr.solvers.linadvec import LinearAdvectionElements


class BaseFluidElements:
    formulations = ['std', 'dual']

    privarmap = {2: ['rho', 'u', 'v', 'p'],
                 3: ['rho', 'u', 'v', 'w', 'p']}

    convarmap = {2: ['rho', 'rhou', 'rhov', 'p'],
                 3: ['rho', 'rhou', 'rhov', 'rhow', 'p']}

    dualcoeffs = convarmap

    visvarmap = {
        2: [('density', ['rho']),
            ('velocity', ['u', 'v']),
            ('pressure', ['p'])],
        3: [('density', ['rho']),
            ('velocity', ['u', 'v', 'w']),
            ('pressure', ['p'])]
    }

    @staticmethod
    def pri_to_con(pris, rhob, cfg):
        rho, p = pris[0], pris[-1]

        # Multiply velocity components by rhob
        rhovs = [rhob*c for c in pris[1:-1]]
        return [rho] + rhovs + [p]

    @staticmethod
    def con_to_pri(cons, cfg):
        ptr = int(len(cons)/2)
        rho, p = cons[0], cons[ptr-1]
        rhob = cons[ptr]

        # Divide momentum components by rhob
        vs = [rhov/rhob for rhov in cons[1:ptr-1]]
        return [rho] + vs + [p]


class LinearEulerElements(BaseFluidElements, LinearAdvectionElements):
    def set_backend(self, *args, **kwargs):
        super().set_backend(*args, **kwargs)

        # Can elide interior flux calculations at p = 0
        if self.basis.order == 0:
            return

        # Register our flux Kernels
        self._be.pointwise.register('pyfr.solvers.lineuler.kernels.tflux')
        self._be.pointwise.register('pyfr.solvers.lineuler.kernels.tfluxlin')

        # Template parameters for the flux Kernels
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'bnvars': self.bnvars,
            'nverts': len(self.basis.linspts),
            'c': self.cfg.items_as('constants', float),
            'jac_exprs': self.basis.jac_exprs
        }

        # Helpers
        c, l = 'curved', 'linear'
        r, s = self._mesh_regions, self._slice_mat


        """I changed it back here"""

        if c in r and 'flux' not in self.antialias:
            self.kernels['tdisf_curved'] = lambda uin: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nupts, r[c]],
                u=s(self.scal_upts[uin], c), f=s(self._vect_upts, c),
                smats=self.curved_smat_at('upts')
            )
        elif c in r:
            self.kernels['tdisf_curved'] = lambda: self._be.kernel(
                'tflux', tplargs=tplargs, dims=[self.nqpts, r[c]],
                u=s(self._scal_qpts, c), f=s(self._vect_qpts, c),
                smats=self.curved_smat_at('qpts')
            )

        if l in r and 'flux' not in self.antialias:
            self.kernels['tdisf_linear'] = lambda uin: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nupts, r[l]],
                u=s(self.scal_upts[uin], l), f=s(self._vect_upts, l),
                verts=self.ploc_at('linspts', l), upts=self.upts
            )
        elif l in r:
            self.kernels['tdisf_linear'] = lambda: self._be.kernel(
                'tfluxlin', tplargs=tplargs, dims=[self.nqpts, r[l]],
                u=s(self._scal_qpts, l), f=s(self._vect_qpts, l),
                verts=self.ploc_at('linspts', l), upts=self.qpts
            )
