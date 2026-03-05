from pathlib import Path

import numpy as np

from pyfr.mpiutil import get_comm_rank_root, mpi
from pyfr.plugins.base import BaseSolverPlugin, RegionMixin, init_csv
from pyfr.util import first


class SFDDampingPlugin(RegionMixin, BaseSolverPlugin):
    name = 'sfd'
    systems = ['*']
    formulations = ['std']
    dimensions = [2, 3]

    def __init__(self, intg, cfgsect):
        super().__init__(intg, cfgsect)

        # Control parameters
        self._chi = self.cfg.getfloat(cfgsect, 'chi')
        self._delta = self.cfg.getfloat(cfgsect, 'delta')

        # Prepare for the error estimation
        self._err_estimate(intg, cfgsect)

        # Solution info
        self._banks = intg.system.ele_banks
        self._rgn_by_doff = {doff: rgn for doff, _, rgn in self._ele_regions}

        # Initialize a base state
        self._init_qbar(intg)

    def _init_qbar(self, intg):
        # Initialization the targeted base flow
        ridx = getattr(intg, '_idxcurr', 0)
        self._qbar = {}
        for doff, _, _ in self._ele_regions:
            self._qbar[doff] = self._banks[doff][ridx].get().copy()

        self._tlast = float(intg.tcurr)
        self._diag_tlast = float(intg.tcurr)

    def _err_estimate(self, intg, cfgsect):
        # MPI info
        comm, rank, root = get_comm_rank_root()

        # Output frequency for error estimation
        self._diag_dt = self.cfg.getfloat(cfgsect, 'diag-dt', 0)
        self._diag_tlast = intg.tcurr

        # Norm used on residual
        self.lp = self.cfg.getfloat(cfgsect, 'norm', 2)

        # Reduction parameters
        if self.lp == float('inf'):
            self._lp_exp = 1
            self._np_op = np.maximum
            self._mpi_op = mpi.MAX
        else:
            self._lp_exp = self.lp
            self._np_op = np.add
            self._mpi_op = mpi.SUM

        intg.call_plugin_dt(intg.tcurr, self._diag_dt)
        if rank == root:
            header = ['t, dt'] + first(intg.system.ele_map.values()).convars
            self.csv = init_csv(self.cfg, cfgsect, ','.join(header), nflush=1)
        else:
            self.csv = None

    def _diag_norms(self, intg):
        if self.csv is None:
            return
        if intg.tcurr - self._diag_tlast < self._diag_dt - self.tol:
            return
        
        # MPI info
        comm, rank, root = get_comm_rank_root()

        # Compute the norms between current and base states
        ridx = getattr(intg, '_idxcurr', 0)

        norms = []
        for doff, rgn in self._rgn_by_doff.items():
            q = self._banks[doff][ridx].get()
            dq = q[..., rgn] - self._qbar[doff][..., rgn]
            dq = dq.swapaxes(0, 1).reshape(self.nvars, -1)
            norms.append(np.linalg.norm(dq, axis=1, ord=self.lp))

        # Reduce over each element type in our domain
        resid = self._np_op.reduce([n**self._lp_exp for n in norms])

        # Reduce over all domains and, if we are the root rank, output
        if rank != root:
            comm.Reduce(resid, None, op=self._mpi_op, root=root)
        else:
            comm.Reduce(mpi.IN_PLACE, resid, op=self._mpi_op, root=root)

            # Post process
            resid = (r**(1 / self._lp_exp) for r in resid)
            
            # Write
            self.csv(intg.tcurr, intg.tcurr - self._tlast, *resid)

        # Update
        self._diag_tlast = float(intg.tcurr)

    def __call__(self, intg):
        dt = float(intg.tcurr - self._tlast)
        if dt <= self.tol:
            return

        # Get scaling parameters
        ridx = getattr(intg, '_idxcurr', 0)
        af = np.exp(-dt / self._delta)
        one_minus_af = 1.0 - af
        ad = np.exp(-self._chi*dt) 

        for doff, rgn in self._rgn_by_doff.items():
            # Get solutions
            q = self._banks[doff][ridx].get()
            qb = self._qbar[doff]

            qv = q[..., rgn]
            qbv = qb[..., rgn]

            # First update the low-pass filtered state using q^n
            qbv *= af
            qbv += one_minus_af*qv

            # Then damp the solution towards qbar^{n+1}.
            qv[:] = qbv + ad*(qv - qbv)

            # update the solution
            self._banks[doff][ridx].set(q)

        # Estimate the error
        intg._invalidate_caches()
        self._diag_norms(intg)

        # Update
        self._tlast = float(intg.tcurr)
