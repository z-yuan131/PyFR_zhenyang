import ctypes as ct
from ctypes.util import find_library
from pathlib import Path
import numpy as np

from pyfr.mpiutil import get_comm_rank_root, mpi
from pyfr.plugins.base import BaseSolverPlugin
# from pyfr.plugins.arpack_test.eigs_parpack import ParpackRCISession

def _sort_idx(vals, which):
    if which in {'LR', 'LA'}:
        key = np.real(vals)
    elif which in {'LM'}:
        key = np.abs(vals)
    elif which in {'LI'}:
        key = np.imag(vals)
    else:
        raise ValueError(f'Unsupported which={which}')
    return np.argsort(key)[::-1]

class _ParpackLib:
    def __init__(self):
        libname = find_library('parpack') or 'libparpack.so'
        self.lib = ct.CDLL(libname)

        self.pdnaupd = self.lib.pdnaupd_c
        self.pdneupd = self.lib.pdneupd_c

        cint = ct.c_int
        cdouble = ct.c_double
        ccharp = ct.c_char_p
        nd_i = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
        nd_d = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
        nd_d_f = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS')

        self.pdnaupd.argtypes = [
            cint,                        # MPI_Fint comm
            ct.POINTER(cint),            # ido
            ccharp,                      # bmat
            cint,                        # n (local row count)
            ccharp,                      # which
            cint,                        # nev
            cdouble,                     # tol
            nd_d,                        # resid
            cint,                        # ncv
            nd_d_f,                      # v
            cint,                        # ldv
            nd_i,                        # iparam
            nd_i,                        # ipntr
            nd_d,                        # workd
            nd_d,                        # workl
            cint,                        # lworkl
            ct.POINTER(cint),            # info
        ]
        self.pdnaupd.restype = None

        self.pdneupd.argtypes = [
            cint,                        # MPI_Fint comm
            cint,                        # rvec
            ccharp,                      # howmny
            nd_i,                        # select
            nd_d,                        # dr
            nd_d,                        # di
            nd_d_f,                      # z
            cint,                        # ldz
            cdouble,                     # sigmar
            cdouble,                     # sigmai
            nd_d,                        # workev
            ccharp,                      # bmat
            cint,                        # n
            ccharp,                      # which
            cint,                        # nev
            cdouble,                     # tol
            nd_d,                        # resid
            cint,                        # ncv
            nd_d_f,                      # v
            cint,                        # ldv
            nd_i,                        # iparam
            nd_i,                        # ipntr
            nd_d,                        # workd
            nd_d,                        # workl
            cint,                        # lworkl
            ct.POINTER(cint),            # info
        ]
        self.pdneupd.restype = None


class ParpackRCISession:
    """Stateful PARPACK reverse-communication session.

    Call `iterate(yloc=...)` once per available operator application.
    When PARPACK needs A*x, `iterate` returns `{'state': 'need_op', 'x': ...}`.
    Provide the propagated vector as the next `yloc`.
    """

    def __init__(self, nloc, nglob, *, k=6, which='LR', ncv=None, tol=0.0,
                 maxiter=None, v0=None, comm=mpi.COMM_WORLD):
        self.nloc = int(nloc)
        self.nglob = int(nglob)
        self.k = int(k)
        self.which = which
        self.comm = comm

        if self.nglob <= 0:
            raise ValueError('Global problem size must be positive')
        if self.k < 1 or self.k >= self.nglob:
            raise ValueError('k must satisfy 1 <= k < n')

        self.ncv = int(ncv or min(max(2*self.k + 1, 20), self.nglob))
        if not (self.k + 1 <= self.ncv <= self.nglob):
            raise ValueError('ncv must satisfy k + 1 <= ncv <= n')

        self.maxiter = int(maxiter or max(300, 20*self.nglob))

        self._parpack = _ParpackLib()

        self._bmat = b'I'
        self._which_b = which.encode()
        self._tol = float(tol)

        self._ido = ct.c_int(0)
        self._info = ct.c_int(0)

        if v0 is not None:
            v0 = np.asarray(v0, dtype=np.float64)
            if len(v0) != self.nloc:
                raise ValueError('v0 must be local vector of length nloc')
            self._info = ct.c_int(1)

        self._resid = np.zeros(self.nloc, dtype=np.float64)
        if v0 is not None:
            self._resid[:] = v0

        self._v = np.zeros((self.nloc, self.ncv), dtype=np.float64, order='F')
        self._iparam = np.zeros(11, dtype=np.int32)
        self._ipntr = np.zeros(14, dtype=np.int32)
        self._workd = np.zeros(3*self.nloc, dtype=np.float64)
        self._lworkl = 3*self.ncv*(self.ncv + 2)
        self._workl = np.zeros(self._lworkl, dtype=np.float64)

        # ARPACK iparam setup (mode 1, exact shifts)
        self._iparam[0] = 1
        self._iparam[2] = self.maxiter
        self._iparam[3] = 1
        self._iparam[6] = 1

        self._commf = ct.c_int(comm.py2f())
        self._ldv = ct.c_int(self.nloc)
        self._n_c = ct.c_int(self.nloc)
        self._nev_c = ct.c_int(self.k)
        self._ncv_c = ct.c_int(self.ncv)
        self._lworkl_c = ct.c_int(self._lworkl)
        self._tol_c = ct.c_double(self._tol)

        self._pending_ido = None
        self._finished = False

    @property
    def done(self):
        return self._finished

    @property
    def ido(self):
        return self._ido.value

    @property
    def nconv(self):
        return int(self._iparam[4])

    def _xptr_yptr(self):
        xptr = self._ipntr[0] - 1
        yptr = self._ipntr[1] - 1
        return xptr, yptr

    def _get_requested_x(self):
        xptr, _ = self._xptr_yptr()
        return self._workd[xptr:xptr + self.nloc].copy()

    def _put_operator_result(self, yloc):
        yloc = np.asarray(yloc, dtype=np.float64)
        if yloc.shape != (self.nloc,):
            raise ValueError('yloc must be local vector of length nloc')

        _, yptr = self._xptr_yptr()
        self._workd[yptr:yptr + self.nloc] = yloc

    def _apply_identity_b(self):
        xptr, yptr = self._xptr_yptr()
        self._workd[yptr:yptr + self.nloc] = self._workd[xptr:xptr + self.nloc]

    def _pdnaupd(self):
        self._parpack.pdnaupd(
            self._commf,
            ct.byref(self._ido),
            self._bmat,
            self._n_c,
            self._which_b,
            self._nev_c,
            self._tol_c,
            self._resid,
            self._ncv_c,
            self._v,
            self._ldv,
            self._iparam,
            self._ipntr,
            self._workd,
            self._workl,
            self._lworkl_c,
            ct.byref(self._info)
        )

        if self._info.value < 0:
            raise RuntimeError(f'pdnaupd_c failed with info={self._info.value}')

    def iterate(self, yloc=None):
        """Advance PARPACK by one reverse-communication exchange.

        Returns:
            {'state': 'need_op', 'ido': ido, 'x': xloc}
            {'state': 'done', 'ido': 99}
        """
        if self._finished:
            return {'state': 'done', 'ido': 99}

        if self._pending_ido in (-1, 1):
            if yloc is None:
                raise ValueError('PARPACK expects operator output yloc for pending ido')
            self._put_operator_result(yloc)
        elif yloc is not None:
            raise ValueError('yloc provided but PARPACK is not waiting for operator output')

        self._pending_ido = None

        while True:
            self._pdnaupd()

            if self._ido.value in (-1, 1):
                self._pending_ido = self._ido.value
                return {
                    'state': 'need_op',
                    'ido': self._ido.value,
                    'x': self._get_requested_x()
                }
            if self._ido.value == 2:
                # bmat='I' so B*x = x; satisfy internally and continue.
                self._apply_identity_b()
                continue
            if self._ido.value == 99:
                self._finished = True
                return {'state': 'done', 'ido': 99}

            raise RuntimeError(f'Unexpected ido={self._ido.value}')

    def extract(self, return_eigenvectors=True):
        if not self._finished:
            raise RuntimeError('Cannot extract eigenpairs before PARPACK reports ido=99')

        nconv = self.nconv
        if nconv <= 0:
            raise RuntimeError('No converged eigenvalues')

        rvec = ct.c_int(1 if return_eigenvectors else 0)
        howmny = b'A'
        select = np.zeros(self.ncv, dtype=np.int32)
        dr = np.zeros(self.k + 1, dtype=np.float64)
        di = np.zeros(self.k + 1, dtype=np.float64)
        z = np.zeros((self.nloc, self.k + 1), dtype=np.float64, order='F')
        workev = np.zeros(3*self.ncv, dtype=np.float64)
        sigmar = ct.c_double(0.0)
        sigmai = ct.c_double(0.0)
        info_e = ct.c_int(0)

        self._parpack.pdneupd(
            self._commf,
            rvec,
            howmny,
            select,
            dr,
            di,
            z,
            ct.c_int(self.nloc),
            sigmar,
            sigmai,
            workev,
            self._bmat,
            self._n_c,
            self._which_b,
            self._nev_c,
            self._tol_c,
            self._resid,
            self._ncv_c,
            self._v,
            self._ldv,
            self._iparam,
            self._ipntr,
            self._workd,
            self._workl,
            self._lworkl_c,
            ct.byref(info_e)
        )

        if info_e.value != 0:
            raise RuntimeError(f'pdneupd_c failed with info={info_e.value}')

        vals = dr[:self.k] + 1j*di[:self.k]
        order = _sort_idx(vals, self.which)
        vals = vals[order]

        if not return_eigenvectors:
            return vals

        vecs_local = z[:, :self.k][:, order]
        return vals, vecs_local




class ArnoldiPlugin(BaseSolverPlugin):
    """PARPACK-based eigen-analysis plugin.

    Workflow:
    1) Every `sample-every` accepted solver steps, feed current state vector to PARPACK.
    2) PARPACK performs one reverse-communication iteration.
    3) Write PARPACK's newest orthogonal vector back to the time-stepper.
    4) Repeat until PARPACK reports ido == 99.
    """

    name = 'arnoldi'
    systems = ['linear-navier-stokes']
    formulations = ['std']
    dimensions = [2, 3]

    def __init__(self, intg, cfgsect):
        super().__init__(intg, cfgsect)

        comm, rank, root = get_comm_rank_root()
        self._comm = comm
        self._rank = rank
        self._root = root

        self._k = self.cfg.getint(cfgsect, 'n-eigs', 6)
        self._which = self.cfg.get(cfgsect, 'which', 'LR')
        self._ncv = self.cfg.getint(cfgsect, 'ncv', 0) or None
        self._tol = self.cfg.getfloat(cfgsect, 'eig-tol', 1e-8)
        self._maxiter = self.cfg.getint(cfgsect, 'maxiter', 0) or None
        self._init_mode = self.cfg.get(cfgsect, 'initial-vector', 'arpack')
        if self._init_mode not in {'arpack', 'current', 'random'}:
            raise ValueError('solver-plugin-arnoldi: initial-vector must be arpack, current, or random')

        # How often to sample the time-stepper and do a PARPACK exchange
        self._sample_dt = self.cfg.getfloat(cfgsect, 'sample-dt')
        self.tout_last = intg.tcurr

        # Register our output times with the integrator
        intg.call_plugin_dt(intg.tcurr, self._sample_dt)

        #self._run_at_step = self.cfg.getint(cfgsect, 'run-at-step', 0)
        self._abort_when_done = self.cfg.getbool(cfgsect, 'abort-when-done', True)

        # The whole naming routine should be checked, making sure they are aligned with the other plugins and when making output for the eigenvectors, we can simply use native writer.
        basedir = Path(self.cfg.getpath(cfgsect, 'basedir', '.', abs=True))
        basename = self.cfg.get(cfgsect, 'basename', 'parpack_lns')

        self._vals_path = Path(self.cfg.getpath(cfgsect, 'eigs-file',
                                str(basedir / f'{basename}.eigs.npy'), abs=True))
        self._vecs_path = Path(self.cfg.getpath(cfgsect, 'vecs-file',
                                str(basedir / f'{basename}.vecs.rank{rank}.npy'), abs=True))

        # Get the element map and register
        self._banks = intg.system.ele_banks
        self._ridx = getattr(intg, '_idxcurr', 0)

        self._parts = []
        for b in self._banks:
            sh = b[self._ridx].ioshape
            n = int(np.prod(sh))
            self._parts.append((sh, n))
        self._nloc = sum(n for _, n in self._parts)
        self._nglob = int(self._comm.allreduce(self._nloc))

        self._session = None
        #self._last_exchange_step = None
        self._done = False

        # If we're not restarting then make sure we initializae the 
        # PARPACK session immediately so that it can start iterating right away
        if not intg.isrestart:
            self.tout_last -= self._sample_dt

    def _pack_reg(self, ridx):
        out = []
        for b in self._banks:
            out.append(b[ridx].get().reshape(-1))

        if not out:
            return np.empty(0, dtype=np.float64)

        return np.concatenate(out).astype(np.float64, copy=False)

    def _unpack_to_reg(self, ridx, xloc):
        xloc = np.asarray(xloc)

        off = 0
        for b, (sh, n) in zip(self._banks, self._parts):
            arr = xloc[off:off + n].reshape(sh)
            b[ridx].set(arr)
            off += n

    def _set_stepper_vector(self, intg, xloc):
        ridx = getattr(intg, '_idxcurr', 0)
        self._unpack_to_reg(ridx, xloc)
        intg._idxcurr = ridx
        intg._invalidate_caches()

    def _save_results(self, vals, vecs_loc):
        self._vals_path.parent.mkdir(parents=True, exist_ok=True)

        if self._rank == self._root:
            np.save(self._vals_path, vals)

        np.save(self._vecs_path, vecs_loc)

    def _global_norm(self, xloc):
        l2loc = float(np.vdot(xloc, xloc).real)
        l2glob = self._comm.allreduce(l2loc)
        return float(np.sqrt(l2glob))

    def _make_initial_vector(self, intg):
        if self._init_mode == 'arpack':
            return None
        if self._init_mode == 'current':
            v0 = self._pack_reg(getattr(intg, '_idxcurr', 0))
        else:
            rng = np.random.default_rng(1234 + self._rank)
            v0 = rng.standard_normal(self._nloc)

        nrm = self._global_norm(v0)

        if nrm == 0:
            raise RuntimeError('Failed to construct a non-zero initial vector for PARPACK')

        return v0 / nrm

    def _start_parpack(self, intg):
        v0 = self._make_initial_vector(intg)

        self._session = ParpackRCISession(
            self._nloc,
            self._nglob,
            k=self._k,
            which=self._which,
            ncv=self._ncv,
            tol=self._tol,
            maxiter=self._maxiter,
            v0=v0,
            comm=self._comm,
        )

        out = self._session.iterate()

        if self._rank == self._root:
            print('Parpack handshake', 'time', intg.tcurr, 'ido', out['ido'], flush=True)


        if out['state'] == 'done':
            vals, vecs_loc = self._session.extract(return_eigenvectors=True)
            self._save_results(vals, vecs_loc)
            self._done = True
            return

        self._set_stepper_vector(intg, out['x'])
        #self._last_exchange_step = intg.nacptsteps

        

    def _parpack_exchange(self, intg):
        yloc = self._pack_reg(getattr(intg, '_idxcurr', 0))
        out = self._session.iterate(yloc=yloc)

        if self._rank == self._root:
            print('Parpack handshake', 'time', intg.tcurr, 'ido', out['ido'], flush=True)

        if out['state'] == 'done':
            vals, vecs_loc = self._session.extract(return_eigenvectors=True)
            self._save_results(vals, vecs_loc)
            self._done = True
            return

        self._set_stepper_vector(intg, out['x'])
        #self._last_exchange_step = intg.nacptsteps

        

    def _parpack_routine(self, intg):
        if self._session is None:
            self._start_parpack(intg)
        else:
            self._parpack_exchange(intg)

    def __call__(self, intg):
        if self._done:
            return

        # I don't think this is necessary since arnoldi algorithm should start immediately, but just in case, we can delay starting until a certain number of accepted steps have passed
        #if intg.nacptsteps < self._run_at_step:
        #    return

        if intg.tcurr - self.tout_last < self._sample_dt - self.tol:
            return

        # Perform the PARPACK handshake
        self._parpack_routine(intg)

        #if self._session is None:
        #    self._start_parpack(intg)
        #elif intg.nacptsteps - self._last_exchange_step >= self._sample_dt:
        #    self._parpack_exchange(intg)

        # I think this should be default, why still running after PARPACK is done? Just abort the time-stepper if PARPACK is done and we want to abort
        if self._done and self._abort_when_done:
            intg.plugin_abort('Arnoldi/PARPACK finished (ido=99)')

        # Update the last output time
        self.tout_last = intg.tcurr
