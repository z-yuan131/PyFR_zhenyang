<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
% for mod, name in src_macros:
    <%include file='${mod}'/>
% endfor

#include <stdio.h>



<%pyfr:kernel name='negdivconf' ndim='2'
              t='scalar fpdtype_t'
              tdivtconf='inout fpdtype_t[${str(nvars)}]'
              ploc='in fpdtype_t[${str(ndims)}]'
              u='in fpdtype_t[${str(nvars)}]'
              rcpdjac='in fpdtype_t'>
fpdtype_t src[${nvars}] = {};

% for mod, name in src_macros:
    ${pyfr.expand(name, 't', 'u', 'ploc', 'src')};
% endfor


% for i in range(nvars):
    tdivtconf[${i}] = -rcpdjac*tdivtconf[${i}] + src[${i}]; // + linsrc[${i}];
% endfor


// printf("%f, %f, %f\n", tdivtconf[2], tdivtconf[1], tdivtconf[0]);

</%pyfr:kernel>
