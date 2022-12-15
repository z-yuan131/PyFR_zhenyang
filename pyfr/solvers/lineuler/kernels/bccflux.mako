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

    // Compute the RHS
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'norm_nl', 'ur')};

    // Perform the Riemann solve
    fpdtype_t fn[${bnvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'norm_nl', 'fn')};

    // Scale and write out the common normal fluxes
% for i in range(bnvars):
    ul[${i}] = mag_nl*fn[${i}];
    // set 0 to all base flow flux
    ul[${i+bnvars}] = 0;
% endfor
</%pyfr:kernel>
