# -*- coding: utf-8 -*-
<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:kernel name='gradcoru' ndim='2'
              gradu='in fpdtype_t[${str(ndims)}][${str(nvars)}]'
              smats='in fpdtype_t[${str(ndims)}][${str(ndims)}]'
              rcpdjac='in fpdtype_t'
              gradbaseu='out fpdtype_t[${str(ndims)}][${str(bnvars)}]'>
    fpdtype_t tmpgradu[][${nvars}] = ${pyfr.array('gradu[{i}][{j}]', i=ndims, j=nvars)};

% for i, j in pyfr.ndrange(ndims, bnvars):
    gradbaseu[${i}][${j}] = rcpdjac*(${' + '.join(f'smats[{k}][{i}]*tmpgradu[{k}][{j+bnvars}]'
                                              for k in range(ndims))});
% endfor
</%pyfr:kernel>
