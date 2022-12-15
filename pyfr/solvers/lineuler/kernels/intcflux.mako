# -*- coding: utf-8 -*-u
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.lineuler.kernels.rsolvers.${rsolver}'/>

<%pyfr:kernel name='intcflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              ur='inout view fpdtype_t[${str(nvars)}]'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};

    // Perform the Riemann solve
    fpdtype_t fn[${bnvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'fn')};

    // Scale and write out the common normal fluxes
% for i in range(bnvars):
    ul[${i}] =  mag_nl*fn[${i}];
    ur[${i}] = -mag_nl*fn[${i}];

    // set base flow flux to 0
    ul[${i+bnvars}] = 0;
    ur[${i+bnvars}] = 0;
% endfor
</%pyfr:kernel>
