<% import numpy as np %>

<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='alm' params='t, u, ploc, src' externs='nloc, forc'>

  fpdtype_t g = 0.0;
  fpdtype_t f[${ndims}] = {};
  fpdtype_t invrho = 1/u[0];

  // Forcing with the Gaussian filter
  % for i in range(npts):

    fpdtype_t r2 = 0.0;
    % for j in range(ndims):
      r2 += (ploc[${j}] - nloc[${i}][${j}])*(ploc[${j}] - nloc[${i}][${j}]);
    % endfor

    g = exp(-r2*${eph2});
     
    % for j in range(ndims):
      f[${j}] += g*forc[${i}][${j}];
    % endfor 

  % endfor

  // Momentum eq.
  % for i in range(ndims):
    src[${i + 1}] += f[${i}] * ${ephpi3};
  % endfor

  // Energy eq.
  % for i in range(ndims):
    src[${nvars - 1}] += f[${i}]*(u[${i+1}] * invrho ) * ${ephpi3};
  % endfor
</%pyfr:macro>
