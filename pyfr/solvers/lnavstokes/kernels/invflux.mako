<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

#include <stdio.h>


<%pyfr:macro name='inviscid_flux' params='s, f' externs='ub'>
    fpdtype_t rhob = ub[0], rho = s[0], invrhob = 1.0/ub[0];
    fpdtype_t pb = ub[${nvars - 1}], p = s[${nvars - 1}];

    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = s[${i + 1}]*rhob + ub[${i + 1}]*rho;
    f[${i}][${nvars - 1}] = ub[${i + 1}]*p + ${c['gamma']}*pb*s[${i + 1}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = s[${i + 1}]*ub[${j + 1}]${' + p*invrhob' if i == j else ''};
% endfor



// printf("%f, %f\n", ub[0], s[0]);


</%pyfr:macro>


