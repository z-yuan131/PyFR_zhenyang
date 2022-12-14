# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.lineuler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.lineuler.kernels.bcs.${bctype}'/>

<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Split perturbation variables
    fpdtype_t ulp[${${bnvars}}]
% for i in range(bnvars):
    ulp[${i}] = ul[${i}];
% endfor

    // Compute the RHS
    fpdtype_t urp[${bnvars}];
    ${pyfr.expand('bc_rsolve_state', 'ulp', 'norm_nl', 'urp')};

    // Perform the Riemann solve
    fpdtype_t fn[${bnvars}];
    ${pyfr.expand('rsolve', 'ulp', 'urp', 'norm_nl', 'fn')};

    // Scale and write out the common normal fluxes
% for i in range(bnvars):
    ul[${i}] = mag_nl*fn[${i}];
% endfor
</%pyfr:kernel>
