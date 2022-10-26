# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
// <%include file='pyfr.solvers.baseadvec.kernels.negdivufb'/>

<%pyfr:kernel name='negdivconflin' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'
              cu='in fpdtype_t[${str(ndims)}]'
              >

    // Compute the C@U term in the formula
    // fpdtype_t ftemp[${nvars}];
    // ${pyfr.expand('utdivtconf', 'u', 'ub', 'divub', 'ftemp')};


    % for i, ex in enumerate(srcex):
        tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + ${ex} - cu[${i}];
    % endfor


</%pyfr:kernel>
