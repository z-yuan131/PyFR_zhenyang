# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.lineuler.kernels.flux'/>

<%pyfr:kernel name='tflux' ndim='2'
              u='in fpdtype_t[${str(nvars)}]'
              f='out fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'>

    // Compute the flux
    fpdtype_t ftemp[${ndims}][${bnvars}];
    fpdtype_t p, v[${ndims}];
    ${pyfr.expand('inviscid_flux', 'u', 'ftemp', 'p', 'v')};

    // Transform the fluxes
% for i, j in pyfr.ndrange(ndims, bnvars):
    f[${i}][${j}] = ${' + '.join(f'smats[{i}][{k}]*ftemp[{k}][{j}]'
                                 for k in range(ndims))};

    // set all base flow flux to 0
    f[${i}][${j+bnvars}] = 0;
% endfor
</%pyfr:kernel>
