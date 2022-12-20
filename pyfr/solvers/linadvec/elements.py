# -*- coding: utf-8 -*-

from pyfr.polys import get_polybasis
from pyfr.solvers.baseadvec import BaseAdvectionElements


class LinearAdvectionElements(BaseAdvectionElements):
    @property
    def _scratch_bufs(self):
        bufs = {'scal_fpts', 'vect_upts', 'vect_fpts','basegrad_upts'}

        if 'flux' in self.antialias:
            bufs |= {'scal_qpts', 'vect_qpts'}


        if self._soln_in_src_exprs:
            bufs |= {'scal_upts_cpy'}

        return bufs

    def set_backend(self, backend, nscalupts, nonce, linoff):
        super().set_backend(backend, nscalupts, nonce, linoff)

        kernels = self.kernels

        # For a linear advection type solver, gradient of baseflow is needed:
        kprefix = 'pyfr.solvers.linadvec.kernels'
        slicem = self._slice_mat


        # Register our pointwise kernels
        self._be.pointwise.register(f'{kprefix}.gradcoru')
        self._be.pointwise.register(f'{kprefix}.gradcorulin')

        # Mesh regions
        regions = self._mesh_regions

        if abs(self.cfg.getfloat('solver-interfaces', 'ldg-beta')) == 0.5:
            kernels['copy_fpts'] = lambda: self._be.kernel(
                'copy', self._vect_fpts.slice(0, self.nfpts), self._scal_fpts
            )

        if self.basis.order > 0:
            kernels['tgradpcoru_upts'] = lambda uin: self._be.kernel(
                'mul', self.opmat('M4 - M6*M0'), self.scal_upts[uin],
                out=self._vect_upts
            )
        kernels['tgradcoru_upts'] = lambda: self._be.kernel(
            'mul', self.opmat('M6'), self._vect_fpts.slice(0, self.nfpts),
            out=self._vect_upts, beta=float(self.basis.order > 0)
        )

        # Template arguments for the physical gradient kernel
        tplargs = {
            'ndims': self.ndims,
            'nvars': self.nvars,
            'bnvars': self.bnvars,
            'nverts': len(self.basis.linspts),
            'jac_exprs': self.basis.jac_exprs
        }

        if 'curved' in regions:
            kernels['gradcoru_upts_curved'] = lambda: self._be.kernel(
                'gradcoru', tplargs=tplargs,
                dims=[self.nupts, regions['curved']],
                gradu=slicem(self._vect_upts, 'curved'),
                smats=self.curved_smat_at('upts'),
                rcpdjac=self.rcpdjac_at('upts', 'curved'),
                gradbaseu=slicem(self._basegrad_upts, 'curved')
            )


        if 'linear' in regions:
            kernels['gradcoru_upts_linear'] = lambda: self._be.kernel(
                'gradcorulin', tplargs=tplargs,
                dims=[self.nupts, regions['linear']],
                gradu=slicem(self._vect_upts, 'linear'),
                upts=self.upts, verts=self.ploc_at('linspts', 'linear'),
                gradbaseu=slicem(self._basegrad_upts, 'linear')
            )

        def gradcoru_fpts():
            nupts, nfpts = self.nupts, self.nfpts
            vupts, vfpts = self._vect_upts, self._vect_fpts

            # Exploit the block-diagonal form of the operator
            muls = [self._be.kernel('mul', self.opmat('M0'),
                           vupts.slice(i*nupts, (i + 1)*nupts),
                           vfpts.slice(i*nfpts, (i + 1)*nfpts))
                    for i in range(self.ndims)]

            return self._be.unordered_meta_kernel(muls)

        kernels['gradcoru_fpts'] = gradcoru_fpts

        # What anti-aliasing options we're running with
        if 'flux' in self.antialias and self.basis.order > 0:
            def gradcoru_qpts():
                nupts, nqpts = self.nupts, self.nqpts
                vupts, vqpts = self._vect_upts, self._vect_qpts

                # Exploit the block-diagonal form of the operator
                muls = [self._be.kernel('mul', self.opmat('M7'),
                                        vupts.slice(i*nupts, (i + 1)*nupts),
                                        vqpts.slice(i*nqpts, (i + 1)*nqpts))
                        for i in range(self.ndims)]

                return self._be.unordered_meta_kernel(muls)

            self.kernels['gradcoru_qpts'] = gradcoru_qpts
