# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='bcconu' ndim='1'
              ulin='in view fpdtype_t[${str(nvars)}]'
              ulout='out view fpdtype_t[${str(nvars)}]'
              nlin='in fpdtype_t[${str(ndims)}]'>
% for i in range(bnvars):
      ulout[${i}] = ulin[${i}];
% endfor
</%pyfr:kernel>
