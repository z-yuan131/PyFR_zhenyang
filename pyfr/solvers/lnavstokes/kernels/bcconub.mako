<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='bc_rsolve_state_bf' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = ul[${i + 1}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:kernel name='bcconub' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              nlin='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nlin[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nlin[{i}]', i=ndims)};

    // Baseflow boundary condition is simply copied from the interior
    ${pyfr.expand('bc_rsolve_state_bf', 'ulin', 'norm_nl', 'ulout')};
</%pyfr:kernel>
