# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='cu' ndim='2'
              t='scalar fpdtype_t'
              ub='in fpdtype_t[${str(nvars)}]'
              u='in fpdtype_t[${str(nvars)}]'
              divub='in fpdtype_t[${str(nvars)}][${str(nvars)}]'
              cu='out fpdtype_t[${str(nvars-1)}]'>

    // Compute the average quantities
    fpdtype_t p;
    fpdtype_t v[${nvars}];
    fpdtype_t invrhob = 1.0/ub[0];

% for i in range(ndims):
    v[${i}] = invrhob*u[${i + 1}];
% endfor


% if ndims == 2:

  cu[0] =  u[0]*(ub[1]*divub[1][1] + ub[2]*divub[1][2]);
  cu[0] += ub[0]*v[1]*(divub[1][1] + divub[2][2]);

  cu[1] =  u[0]*(ub[1]*divub[2][1] + ub[2]*divub[2][2]);
  cu[1] += ub[0]*v[2]*(divub[1][1] + divub[2][2]);

  cu[2] =  (${c['gamma'] - 1})*u[3]*(divub[1][1] + divub[2][2]);
  cu[2] += (${1 - c['gamma']})*(v[1]*divub[3][1] + v[2]*divub[3][2]);


% elif ndims == 3:

  cu[0] =  u[0]*(ub[1]*divub[1][1] + ub[2]*divub[1][2] + ub[3]*divub[1][3]);
  cu[0] += ub[0]*v[1]*(divub[1][1] + divub[2][2] + divub[3][3]);

  cu[1] =  u[0]*(ub[1]*divub[2][1] + ub[2]*divub[2][2] + ub[3]*divub[2][3]);
  cu[1] += ub[0]*v[2]*(divub[1][1] + divub[2][2] + divub[3][3]);

  cu[2] =  u[0]*(ub[1]*divub[3][1] + ub[2]*divub[3][2] + ub[3]*divub[3][3]);
  cu[2] += ub[0]*v[3]*(divub[1][1] + divub[2][2] + divub[3][3]);

  cu[3] =  (${c['gamma'] - 1})*u[4]*(divub[1][1] + divub[2][2] + divub[3][3]);
  cu[3] += (${1 - c['gamma']})*(v[1]*divub[4][1] + v[2]*divub[4][2] + v[3]*divub[4][3]);

% endif



</%pyfr:kernel>
