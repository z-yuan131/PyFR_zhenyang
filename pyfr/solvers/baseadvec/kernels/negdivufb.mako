# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='utdivtconf' params='s, sb, divsb ,f'>
    // Compute the average quantities
    fpdtype_t vb[${nvars}];
    fpdtype_t invrhob = 1.0/sb[0];

% for i in range(ndims):
    vb[${i}] = invrhob*sb[${i + 1}];
% endfor

    f[0] = 0;
    f[${nvars - 1}] = 0;
% for i,j in pyfr.ndrange(ndims):
    f[${i+1}] += s[0]*sb[${j+1}]*divsb[${i}][${j}] + s[${j+1}]*divsb[${i}][${j}];
% endfor


</%pyfr:macro>
