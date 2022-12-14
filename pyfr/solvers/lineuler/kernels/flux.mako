# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='u, f, p, v'>

% for i in range(ndims):
    v[${i}] = u[${i+1}]/u[${bnvars}];
% endfor

  // Density and pressure flux
% for i in range(ndims):
    f[${i}][0] = u[0]*u[${i+1+bnvars}] + u[${i+1}];
    f[${i}][${bnvars - 1}] = ${c['gamma']}*u[${nvars - 1}]*v[${i}] + u[${i+1+bnvars}]*u[${bnvars - 1}];
% endfor

  // Notation for the pressure
  p = u[${bnvars - 1}];

  // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = u[${j+1}]*u[${i+1+bnvars}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
