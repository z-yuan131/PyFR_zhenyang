import numpy as np

from pyfr.cache import memoize
from pyfr.readers.native import NativeReader
from pyfr.solvers.lnavstokes.elements import LinearNavierStokesElements
#from pyfr.solvers.navstokes.elements import NavierStokesElements
from pyfr.solvers.lnavstokes.inters import (LinearNavierStokesBaseBCInters,
                                            LinearNavierStokesIntInters,
                                            LinearNavierStokesMPIInters)
#from pyfr.solvers.navstokes.inters import (NavierStokesBaseBCInters,
#                                            NavierStokesIntInters,
#                                            NavierStokesMPIInters)
#from pyfr.solvers.navstokes.system import NavierStokesSystem
from pyfr.solvers.baseadvecdiff import BaseAdvectionDiffusionSystem



class LinearNavierStokesSystem(BaseAdvectionDiffusionSystem):
    name = 'linear-navier-stokes'
    elementscls = LinearNavierStokesElements
    intinterscls = LinearNavierStokesIntInters
    mpiinterscls = LinearNavierStokesMPIInters
    bbcinterscls = LinearNavierStokesBaseBCInters
    #elementscls = NavierStokesElements
    #intinterscls = NavierStokesIntInters
    #mpiinterscls = NavierStokesMPIInters
    #bbcinterscls = NavierStokesBaseBCInters

    def __init__(self, backend, mesh, initsoln, nregs, cfg, serialiser):
        self._mesh = mesh
        self._initsoln = initsoln
        self._m0 = {}

        self.baseflow_upts = []
        self.baseflow_fpts = []
        self.baseflow_grad_upts = []
        self.baseflow_grad_fpts = []

        super().__init__(backend, mesh, initsoln, nregs, cfg, serialiser)

    def commit(self):
        # Build kernels once; they bind to placeholder matrices.
        super().commit()

        # Fill placeholders using this system pipeline.
        self._init_baseflow_cache()

    def _load_eles(self, mesh, initsoln, nregs, nonce):
        eles, elemap = super()._load_eles(mesh, initsoln, nregs, nonce)

        # Cache interpolation operators needed after BaseSystem.commit()
        # deletes self.ele_map.
        self._m0 = {etype: ele.basis.m0 for etype, ele in elemap.items()}

        self._alloc_baseflow_placeholders(elemap, nonce)
        return eles, elemap

    def _alloc_baseflow_placeholders(self, elemap, nonce):
        ndims = self.ndims
        for etype, ele in elemap.items():
            nupts, nvars, neles = ele.nupts, ele.nvars, ele.neles
            nfpts = ele.nfpts

            zub_u = np.zeros((nupts, nvars, neles), dtype=self.backend.fpdtype)
            zub_f = np.zeros((nfpts, nvars, neles), dtype=self.backend.fpdtype)
            zgb_u = np.zeros((ndims, nupts, nvars, neles),
                             dtype=self.backend.fpdtype)
            zgb_f = np.zeros((ndims, nfpts, nvars, neles),
                             dtype=self.backend.fpdtype)

            ub_u = self.backend.matrix((nupts, nvars, neles),
                                       initval=zub_u, tags={'align'},
                                       extent=nonce + 'baseflow_upts')
            ub_f = self.backend.matrix((nfpts, nvars, neles),
                                       initval=zub_f, tags={'align'},
                                       extent=nonce + 'baseflow_fpts')
            gb_u = self.backend.matrix((self.ndims, nupts, nvars, neles),
                                       initval=zgb_u, tags={'align'},
                                       extent=nonce + 'baseflow_grad_upts')
            gb_f = self.backend.matrix((self.ndims, nfpts, nvars, neles),
                                       initval=zgb_f, tags={'align'},
                                       extent=nonce + 'baseflow_grad_fpts')

            self.baseflow_upts.append(ub_u)
            self.baseflow_fpts.append(ub_f)
            self.baseflow_grad_upts.append(gb_u)
            self.baseflow_grad_fpts.append(gb_f)

            # Expose baseflow storage on elements for inter getter methods.
            ele._baseflow_upts = ub_u
            ele._baseflow_fpts = ub_f
            ele._baseflow_grad_upts = gb_u
            ele._baseflow_grad_fpts = gb_f

            # Setup external variables for the source terms in LNS 
            ele._set_external(
                'ub', f'in fpdtype_t[{nvars}]', value=ub_u
            )
            ele._set_external(
                'gb', f'in fpdtype_t[{ndims}][{nvars}]', value=gb_u
            )

    @memoize
    def _compute_grads_graph_baseflow(self, uinbank):
        m = self._mpireqs
        k, *_ = self._get_kernels(uinbank, None)

        def deps(dk, *names):
            return self._kdeps(k, dk, *names)

        g1 = self.backend.graph()
        g1.add_mpi_reqs(m['scal_fpts_recv'])

        # Interpolate solution to flux points
        g1.add_all(k['eles/disu'])

        # Pack/send to MPI neighbors
        g1.add_all(k['mpiint/scal_fpts_pack'], deps=k['eles/disu'])
        for send, pack in zip(m['scal_fpts_send'], k['mpiint/scal_fpts_pack']):
            g1.add_mpi_req(send, deps=[pack])

        for l in k['eles/copy_fpts']:
            g1.add(l, deps=deps(l, 'eles/disu'))
        kdeps = k['eles/copy_fpts'] or k['eles/disu']
        g1.add_all(k['iint/con_u'], deps=kdeps)

        # IMPORTANT: skip BC con_u for baseflow debug path
        g1.add_all(k['bcint/con_ub'], deps=kdeps)

        # Partial corrected gradient
        g1.add_all(k['eles/tgradpcoru_upts'])
        g1.commit()

        g2 = self.backend.graph()

        # MPI con_u
        g2.add_all(k['mpiint/scal_fpts_unpack'])
        for l in k['mpiint/con_u']:
            g2.add(l, deps=deps(l, 'mpiint/scal_fpts_unpack'))

        # Corrected gradient + physical gradient
        g2.add_all(k['eles/tgradcoru_upts'], deps=k['mpiint/con_u'])
        for l in k['eles/gradcoru_u']:
            g2.add(l, deps=deps(l, 'eles/tgradcoru_upts'))

        g2.commit()

        return g1, g2

    def compute_grads_baseflow(self, t, uinbank):
        self._prepare_kernels(t, uinbank, None)
        for graph in self._compute_grads_graph_baseflow(uinbank):
            graph.run()

    def _init_baseflow_cache(self):
        baseflow_soln = self._get_baseflow_soln()

        # Use scratch register if available.
        breg = 1 if self.nregs > 1 else 0
        saved = None
        if breg == 0:
            saved = self.ele_scal_upts(0)

        # Load baseflow into selected register.
        for etype, ebank in zip(self.ele_types, self.ele_banks):
            ebank[breg].set(baseflow_soln[etype])

        # Corrected gradients from the FR pipeline.
        self.preproc(0.0, breg)
        self.compute_grads_baseflow(0.0, breg)

        # Pull baseflow and corrected conservative gradients.
        bf_cons = self.ele_scal_upts(breg)
        grad_cons = [m.get() for m in self.eles_vect_upts]

        # Convert corrected gradients to primitive form for LNS kernels.
        baseflow_upts, basegrad_upts = [], []
        baseflow_fpts, basegrad_fpts = [], []
        for bf, bg in zip(bf_cons, grad_cons):
            _bf = bf.swapaxes(0, 1)
            _bg = np.rollaxis(bg, 2)

            pgrad = self.elementscls.grad_con_to_pri_ns(_bf, _bg, self.cfg)
            # (nvars, ndims, nupts, neles) -> (ndims, nupts, nvars, neles)
            pgrad = np.transpose(np.array(pgrad), (1, 2, 0, 3))

            _bf = self.elementscls.con_to_pri(_bf, self.cfg)
            # (nvars, nupts, neles) -> (nupts, nvars, neles)
            _bf = np.transpose(np.array(_bf), (1, 0, 2))

            baseflow_upts.append(_bf)
            basegrad_upts.append(pgrad)

        # Interpolate cached upts baseflow data to fpts for interfaces.
        for etype, bf, bg in zip(self.ele_types, baseflow_upts, basegrad_upts):
            m0 = self._m0[etype]
            baseflow_fpts.append(np.einsum('ij, jkl -> ikl', m0, bf)) 
            basegrad_fpts.append(np.einsum('ij, ljkm -> likm', m0, bg)) 

        # Update placeholders in-place.
        for ub, s in zip(self.baseflow_upts, baseflow_upts):
            ub.set(s)
        for ub, s in zip(self.baseflow_fpts, baseflow_fpts):
            ub.set(s)
        for gb, g in zip(self.baseflow_grad_upts, basegrad_upts):
            gb.set(g)
        for gb, g in zip(self.baseflow_grad_fpts, basegrad_fpts):
            gb.set(g)

        self.backend.wait()

        if saved is not None:
            for ebank, s in zip(self.ele_banks, saved):
                ebank[0].set(s)

        # Optional debug plotting hook:
        # self._check_baseflow(baseflow_upts, basegrad_upts)

    def _check_baseflow(self, baseflow, basegrad):
        # Do a plot of the baseflow and its gradients to sanity check things before we proceed with the LNS implementation.
        import matplotlib.pyplot as plt
        import matplotlib.tri as tri
        ploc = self.ele_ploc_upts

        mesh, soln, bgrad = [], [], []
        for mm, ss, bg in zip(ploc, baseflow, basegrad):
            m, s = mm.swapaxes(0, 1), ss.swapaxes(0, 1)
            gg = bg.swapaxes(1, 2)
            #m, s, gg = mm, ss, bg
            print(mm.shape, ss.shape, bg.shape)
            print(np.max(gg), np.min(gg))

            mesh.append(m.reshape(2, -1))
            soln.append(s.reshape(4, -1))
            bgrad.append(gg.reshape(2, 4, -1))
            #plt.scatter(m[0], m[1], c=s[1], s=20, edgecolors='none')

        mesh = np.concatenate(mesh, axis=1)
        soln = np.concatenate(soln, axis=1)
        bgrad = np.concatenate(bgrad, axis=2)

        # unique nodes only for plotting
        m = np.unique(mesh, axis=1)
        s = soln[1, np.unique(mesh, axis=1, return_index=True)[1]]
        bg = bgrad[:, :, np.unique(mesh, axis=1, return_index=True)[1]]
        print(bg.shape)
        bg = bg[0, -1]  # 2,4,neles

        plt.figure()
        triang = tri.Triangulation(m[0], m[1])
        plt.tricontourf(triang, s, levels=50)
        plt.colorbar()
        plt.savefig('baseflow.png', dpi=300)

        plt.figure()
        plt.tricontourf(triang, bg, levels=50, cmap = 'bwr')
        plt.colorbar()
        plt.savefig('baseflow_grad.png', dpi=300)
        raise RuntimeError 

    def _get_baseflow_soln(self):
        """Get baseflow solution from restart input or a configured file."""

        ssect = 'solver'
        bsoln = self.cfg.getpath(ssect, 'baseflow-soln', abs=True)
        #bpname = self.cfg.get(ssect, 'baseflow-pname', None)

        #reader = NativeReader(self._mesh.fname, pname=bpname)
        reader = NativeReader(self._mesh.fname)
        try:
            return reader.load_soln(bsoln)
        finally:
            reader.close()


        """There should a section dealing with the situation where baseflow is provided by the configuration expression"""
    


        if not self.cfg.hasopt(ssect, 'baseflow-soln'):
            raise RuntimeError('No baseflow provided; set restart solution or '
                               '[solver-linear-navier-stokes] baseflow-soln')
