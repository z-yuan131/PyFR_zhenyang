# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='u, ub, f, p, v'>

  % for i in range(ndims):
    v[${i}] = u[${i+1}]/ub[0];
  % endfor

  // Density and energy flux
  % for i in range(ndims):
      f[${i}][0] = u[0]*ub[${i+1}] + u[${i+1}];
      f[${i}][${nvars - 1}] = ${c['gamma']}*ub[${nvars - 1}]*v[${i}] + ub[${i+1}]*u[${nvars - 1}];
  % endfor

  // Notation for the pressure
  p = u[${nvars - 1}];

  // Momentum fluxes
  % for i, j in pyfr.ndrange(ndims, ndims):
      f[${i}][${j + 1}] = u[${j+1}]*ub[${i+1}]${' + p' if i == j else ''};
  % endfor

</%pyfr:macro>
