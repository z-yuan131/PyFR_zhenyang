import ctypes as ct
import os
from pathlib import Path
import numpy as np

from pyfr.ctypesutil import LibWrapper
from pyfr.inifile import Inifile
from pyfr.mpiutil import get_comm_rank_root, mpi
from pyfr.plugins.base import BaseSolverPlugin
from pyfr.util import first
from pyfr.writers.native import NativeWriter
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

class _ParpackLib(LibWrapper):
    _libname = 'parpack'

    cint = ct.c_int
    cdouble = ct.c_double
    ccharp = ct.c_char_p
    nd_i = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS')
    nd_d = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    nd_d_f = np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS')

    _functions = [
        (None, 'pdnaupd_c',
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
         ct.POINTER(cint)),           # info
        (None, 'pdneupd_c',
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
         ct.POINTER(cint))            # info
    ]


class ParpackRCISession:
    """Stateful PARPACK reverse-communication session.

    Call `iterate(yloc=...)` once per available operator application.
    When PARPACK needs A*x, `iterate` returns `{'state': 'need_op', 'x': ...}`.
    Provide the propagated vector as the next `yloc`.
    """

    def __init__(self, nloc, nglob, *, k=6, which='LM', ncv=None, tol=1e-6,
                 maxiter=None, v0=None, comm=mpi.COMM_WORLD):
        self.nloc = int(nloc)
        self.nglob = int(nglob)
        self.k = int(k)
        self.which = which
        self.comm = comm

        if self.nglob <= 0:
            # Not meaningful 
            raise ValueError('Global problem size must be positive')
        if self.k < 1 or self.k >= self.nglob:
            raise ValueError('k must satisfy 1 <= k < n')

        # Krylov subspace dimension
        self.ncv = int(ncv or max(2*self.k + 10, 20))
        if not (self.k + 1 <= self.ncv <= self.nglob):
            raise ValueError('ncv must satisfy k + 1 <= ncv <= n')

        self.maxiter = int(maxiter or max(300, 40*self.ncv))

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
    def needs_operator(self):
        return self._pending_ido in (-1, 1)

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
        self._parpack.pdnaupd_c(
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

        self._parpack.pdneupd_c(
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

    def save_state(self, path, *, tcurr=None, nacptsteps=None):
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True) #This will raise error

        payload = {
            'version': np.array(1, dtype=np.int32),
            'nloc': np.array(self.nloc, dtype=np.int64),
            'nglob': np.array(self.nglob, dtype=np.int64),
            'k': np.array(self.k, dtype=np.int64),
            'which': np.array(self.which, dtype='S'),
            'ncv': np.array(self.ncv, dtype=np.int64),
            'tol': np.array(self._tol, dtype=np.float64),
            'maxiter': np.array(self.maxiter, dtype=np.int64),
            'ido': np.array(self._ido.value, dtype=np.int32),
            'info': np.array(self._info.value, dtype=np.int32),
            'pending_ido': np.array(-99 if self._pending_ido is None
                                    else self._pending_ido, dtype=np.int32),
            'finished': np.array(int(self._finished), dtype=np.int8),
            'resid': self._resid,
            'v': self._v,
            'iparam': self._iparam,
            'ipntr': self._ipntr,
            'workd': self._workd,
            'workl': self._workl
        }
        if tcurr is not None:
            payload['tcurr'] = np.array(float(tcurr), dtype=np.float64)
        if nacptsteps is not None:
            payload['nacptsteps'] = np.array(int(nacptsteps), dtype=np.int64)

        tmppath = path.with_suffix(path.suffix + '.tmp')
        with tmppath.open('wb') as f:
            np.savez(f, **payload)
        tmppath.replace(path)

    @classmethod
    def from_state(cls, path, *, comm=mpi.COMM_WORLD):
        with np.load(path, allow_pickle=False) as data:
            if int(data['version']) != 1:
                raise RuntimeError('Unsupported PARPACK checkpoint version')

            sess = cls(
                int(data['nloc']),
                int(data['nglob']),
                k=int(data['k']),
                which=data['which'].tobytes().decode(),
                ncv=int(data['ncv']),
                tol=float(data['tol']),
                maxiter=int(data['maxiter']),
                v0=None,
                comm=comm
            )

            sess._ido.value = int(data['ido'])
            sess._info.value = int(data['info'])

            pending = int(data['pending_ido'])
            sess._pending_ido = None if pending == -99 else pending
            sess._finished = bool(int(data['finished']))

            sess._resid[:] = data['resid']
            sess._v[:] = data['v']
            sess._iparam[:] = data['iparam']
            sess._ipntr[:] = data['ipntr']
            sess._workd[:] = data['workd']
            sess._workl[:] = data['workl']

        return sess




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

        self._k = self.cfg.getint(cfgsect, 'n-eigs')
        self._which = self.cfg.get(cfgsect, 'which', 'LM')
        self._ncv = self.cfg.getint(cfgsect, 'ncv', 0) or None
        self._tol = self.cfg.getfloat(cfgsect, 'eig-tol', 1e-6)
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
        
        # Why we need this? Shouldn't we alway abort when done??
        #self._abort_when_done = self.cfg.getbool(cfgsect, 'abort-when-done', True)

        # The whole naming routine should be checked, making sure they are aligned with the other plugins and when making output for the eigenvectors, we can simply use native writer.
        basedir = self.cfg.getpath(cfgsect, 'basedir', '.', abs=True)
        basename = self.cfg.get(cfgsect, 'basename')

        #self._vals_path = Path(self.cfg.getpath(cfgsect, 'eigs-file',
        #                        str(basedir / f'{basename}.eigs.npy'), abs=True))
        
        self._vals_path = os.path.join(basedir, f'{basename}_eigs.npy')

        # Always write the eigenvectors 
        # self._write_vecs_pyfrs = self.cfg.getbool(cfgsect, 'write-vecs-pyfrs', True)
        #vecs_pyfrs = Path(self.cfg.getpath(cfgsect, 'vecs-pyfrs-file',
        #                   str(basedir / f'{basename}.vecs.pyfrs'), abs=True))
        
        # If restarting from a krylov subspace
        self._resume = self.cfg.getbool(cfgsect, 'resume', False)
        self._state_write_every = self.cfg.getint(cfgsect,
                                                  'state-write-every', 1)
        if self._state_write_every < 1:
            raise ValueError('solver-plugin-arnoldi: state-write-every must '
                             'be >= 1')
        #self._state_path = Path(self.cfg.getpath(
        #    cfgsect, 'state-file',
        #    str(basedir / f'{basename}.state.rank{rank}.npz'), abs=True
        #))
        self._state_path = os.path.join(basedir, f'{basename}.state.rank{rank}.npz')

        # Get the element map and register
        self._ele_types = list(intg.system.ele_types)
        self._emap = intg.system.ele_map
        self._banks = intg.system.ele_banks
        self._ridx = getattr(intg, '_idxcurr', 0)
        self._convars = list(first(self._emap.values()).convars)

        # Get local and global DOF
        self._parts = []
        for etype, b in zip(self._ele_types, self._banks):
            sh = b[self._ridx].ioshape
            n = int(np.prod(sh))
            self._parts.append((etype, sh, n))
        self._nloc = sum(n for _, _, n in self._parts)
        self._nglob = int(self._comm.allreduce(self._nloc))

        # Prepare a writer for the eigenvalues and eigenvectors
        self._vecs_writer = None
        ershapes, erdata = {}, {}
        for etype in self._ele_types:
            if etype in intg.system.mesh.eidxs:
                ershapes[etype] = (self.nvars*self._k, self._emap[etype].nupts)
                erdata[etype] = intg.system.mesh.eidxs[etype]

        self._vecs_writer = NativeWriter.from_integrator(intg, basedir, basename, 
                                                        'arnoldi')
        self._vecs_writer.set_shapes_eidxs(ershapes, erdata)
        self._vecs_async_timeout = self.cfg.getfloat(cfgsect, 'async-timeout', 0)
        """
        if self._write_vecs_pyfrs:
            ershapes, erdata = {}, {}
            for etype in self._ele_types:
                if etype in intg.system.mesh.eidxs:
                    ershapes[etype] = (self.nvars*self._k, self._emap[etype].nupts)
                    erdata[etype] = intg.system.mesh.eidxs[etype]

            self._vecs_writer = NativeWriter.from_integrator(
                intg, vecs_pyfrs.parent, vecs_pyfrs.name, 'arnoldi'
            )
            self._vecs_writer.set_shapes_eidxs(ershapes, erdata)
            self._vecs_async_timeout = self.cfg.getfloat(cfgsect, 'async-timeout', 0)
        """

        self._session = None
        #self._last_exchange_step = None
        self._done = False
        self._last_state_write_step = intg.nacptsteps

        # If restarting, restore the full PARPACK state and requested x.
        if intg.isrestart and self._resume and os.path.exists(self._state_path):
            self._restore_state(intg)
        elif not intg.isrestart:
            # If we're not restarting then make sure we initializae the
            # PARPACK session immediately so that it can start iterating right away.
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
        for b, (_, sh, n) in zip(self._banks, self._parts):
            arr = xloc[off:off + n].reshape(sh)
            b[ridx].set(arr)
            off += n

    def _set_stepper_vector(self, intg, xloc):
        ridx = getattr(intg, '_idxcurr', 0)
        self._unpack_to_reg(ridx, xloc)
        intg._idxcurr = ridx
        intg._invalidate_caches()

    def _prepare_vecs_metadata(self, intg):
        comm, rank, root = get_comm_rank_root()

        fields = [f'ev{j}_{v}' for j in range(self._k) for v in self._convars]

        stats = Inifile()
        stats.set('data', 'fields', ','.join(fields))
        stats.set('data', 'prefix', 'arnoldi')
        stats.set('arnoldi', 'n-eigs', self._k)
        stats.set('arnoldi', 'which', self._which)
        intg.collect_stats(stats)

        if rank == root:
            metadata = {**intg.cfgmeta, 'stats': stats.tostr(),
                        'mesh-uuid': intg.mesh_uuid}
        else:
            metadata = None

        sdata = intg.serialiser.serialise()
        if rank == root:
            metadata |= sdata

        return metadata

    def _prepare_vecs_data(self, vecs_loc):
        data = {}

        # For each element type, stack [mode0(all vars), mode1(all vars), ...]
        # into writer shape (neles, nvars*k, nupts).
        off = 0
        for etype, sh, n in self._parts:
            modes = []
            for j in range(self._k):
                vj = vecs_loc[off:off + n, j].reshape(sh)
                modes.append(vj.transpose(2, 1, 0))

            data[etype] = np.concatenate(modes, axis=1)
            off += n

        return data

    def _save_results(self, intg, vals, vecs_loc):
        # self._vals_path.parent.mkdir(parents=True, exist_ok=True)

        # Save the eigenvalues 
        if self._rank == self._root:
            np.save(self._vals_path, vals)

        #if self._save_vecs_npy:
        #    np.save(self._vecs_path, vecs_loc)

        #if self._vecs_writer is not None:
        #    metadata = self._prepare_vecs_metadata(intg)
        #    data = self._prepare_vecs_data(vecs_loc)
        #    self._vecs_writer.write(data, intg.tcurr, metadata,
        #                            self._vecs_async_timeout)

        # Save eigenvectors using the pyfr native writer
        metadata = self._prepare_vecs_metadata(intg)
        data = self._prepare_vecs_data(vecs_loc)
        self._vecs_writer.write(data, intg.tcurr, metadata,
                                self._vecs_async_timeout)

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
            self._save_results(intg, vals, vecs_loc)
            self._done = True
            return

        self._set_stepper_vector(intg, out['x'])
        #self._last_exchange_step = intg.nacptsteps

    def _save_state(self, intg):
        if self._session is None:
            return

        self._session.save_state(self._state_path, tcurr=intg.tcurr,
                                 nacptsteps=intg.nacptsteps)
        self._last_state_write_step = intg.nacptsteps

    def _restore_state(self, intg):
        self._session = ParpackRCISession.from_state(self._state_path,
                                                     comm=self._comm)
        self._done = self._session.done

        if self._session.done:
            return

        # Sanity check restart time against checkpoint time.
        # This avoids silently resuming with a mismatched solver snapshot.
        # Only root checks and broadcasts the decision.
        chk_ok = True
        msg = None
        if self._rank == self._root:
            with np.load(self._state_path, allow_pickle=False) as data:
                if 'tcurr' in data:
                    tchk = float(data['tcurr'])
                    if abs(tchk - float(intg.tcurr)) > max(self.tol, 1.0e-13):
                        chk_ok = False
                        msg = (f'Arnoldi restart mismatch: checkpoint tcurr={tchk} '
                               f'but solution tcurr={float(intg.tcurr)}')
        chk_ok = self._comm.bcast(chk_ok, root=self._root)
        msg = self._comm.bcast(msg, root=self._root)
        if not chk_ok:
            raise RuntimeError(msg)

        # If PARPACK is waiting for y = A*x, put x back in the stepper.
        if self._session.needs_operator:
            self._set_stepper_vector(intg, self._session._get_requested_x())
        else:
            out = self._session.iterate()
            if out['state'] == 'done':
                self._done = True
            else:
                self._set_stepper_vector(intg, out['x'])

    def _parpack_exchange(self, intg):
        yloc = self._pack_reg(getattr(intg, '_idxcurr', 0))
        out = self._session.iterate(yloc=yloc)

        if self._rank == self._root:
            print('Parpack handshake', 'time', intg.tcurr, 'ido', out['ido'], flush=True)

        if out['state'] == 'done':
            vals, vecs_loc = self._session.extract(return_eigenvectors=True)
            self._save_results(intg, vals, vecs_loc)
            self._done = True
            return

        self._set_stepper_vector(intg, out['x'])
        #self._last_exchange_step = intg.nacptsteps

        

    def _parpack_routine(self, intg):
        if self._session is None:
            self._start_parpack(intg)
        elif self._session.needs_operator:
            self._parpack_exchange(intg)
        else:
            out = self._session.iterate()
            if self._rank == self._root:
                print('Parpack handshake', 'time', intg.tcurr,
                      'ido', out['ido'], flush=True)

            if out['state'] == 'done':
                vals, vecs_loc = self._session.extract(return_eigenvectors=True)
                self._save_results(intg, vals, vecs_loc)
                self._done = True
                return

            self._set_stepper_vector(intg, out['x'])

    def __call__(self, intg):
        if self._done:
            return

        # For now let's don't consider the case where we need to save states
        #if intg.nacptsteps - self._last_state_write_step >= self._state_write_every:
        #    self._save_state(intg)

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

    def finalise(self, intg):
        # self._save_state(intg) # why saves the state at the end?

        super().finalise(intg)

        if self._vecs_writer is not None:
            self._vecs_writer.flush()
