# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='cu' ndim='2'
              t='scalar fpdtype_t'
              ub='in fpdtype_t[${str(nvars)}]'
              u='in fpdtype_t[${str(nvars)}]'
              divub='in fpdtype_t[${str(nvars)}][${str(nvars)}]'
              cu='out fpdtype_t[${str(nvars)}]'>

    // Compute the average quantities
    fpdtype_t vb[${nvars}];
    fpdtype_t v[${nvars}];
    fpdtype_t invrhob = 1.0/ub[0];

% for i in range(ndims):
    v[${i}] = invrhob*u[${i + 1}];
    vb[${i}] = invrhob*ub[${i + 1}];
% endfor


% for i in range(nvars):
    cu[${i+1}] = 1
% endfor


</%pyfr:kernel>
