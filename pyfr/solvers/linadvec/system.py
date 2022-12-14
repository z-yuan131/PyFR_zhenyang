# -*- coding: utf-8 -*-

from pyfr.solvers.baseadvec import BaseAdvectionSystem
from pyfr.util import memoize


class LinearAdvectionSystem(BaseAdvectionSystem):
    @memoize
    def _base_grads_graph(self, uinbank):
        m = self._mpireqs
        k, _ = self._get_kernels(uinbank, None)

        def deps(dk, *names): return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # Interpolate the solution to the flux points
        g1.add_all(k['eles/disub'])

        # Pack and send these interpolated solutions to our neighbours
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disub'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        # Compute the common solution at our internal/boundary interfaces
        for l in k['eles/copy_fpts']:
            g1.add(l, deps=deps(l, 'eles/disub'))
        kdeps = k['eles/copy_fpts'] or k['eles/disub']
        g1.add_all(k['iint/con_u'], deps=kdeps)
        g1.add_all(k['bcint/con_u'], deps=kdeps)

        # Compute the transformed gradient of the partially corrected solution
        g1.add_all(k['eles/tgradpcoru_upts'])
        g1.commit()

        g2 = self.backend.graph()

        # Compute the common solution at our MPI interfaces
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/con_u']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Compute the transformed gradient of the corrected solution
        g2.add_all(k['eles/tgradcoru_upts'], deps=k['mpiint/con_u'])



        # Obtain the physical gradients at the solution points
        for l in k['eles/gradcoru_upts_curved']:
            g2.add(l, deps=deps(l, 'eles/tgradcoru_upts'))
        for l in k['eles/gradcoru_upts_linear']:
            g2.add(l, deps=deps(l, 'eles/tgradcoru_upts'))
        g2.commit()
        #print('_base_grads_graph')
        return g1, g2
