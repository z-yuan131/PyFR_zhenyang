# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<% gamma = c['gamma'] %>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>

    fpdtype_t V_n0 = ${' + '.join('ul[{1}]*nl[{0}]'.format(i + bnvars, i + 1)
                                 for i in range(ndims))};

    fpdtype_t V_n = ${' + '.join('ul[{1}]*nl[{0}]'.format(i, i + 1)
                                 for i in range(ndims))};


    fpdtype_t inv = 1.0 / ul[${bnvars}];
    fpdtype_t c = sqrt(${gamma}*ul[${nvars-1}]*inv);

    fpdtype_t h1 = (V_n0 > 0)
                 ? ul[0] - ul[${bnvars-1}] / (c * c)
                 : 0;

    fpdtype_t h4 = (V_n0 > c)
                 ? V_n/2.0 - ul[${bnvars-1}]/(2.0*c)
                 : 0;

    fpdtype_t h5 = (V_n0 + c >0)
                 ? V_n/2.0 + ul[${bnvars-1}]/(2.0*c)
                 : 0;

    ur[0] = h1 + (h5 - h4)/c;

    ur[${bnvars-1}] = c * (h5 - h4);

% for i in range(ndims):
    ur[${i + 1}] = (h4 + h5)*nl[${i}];
% endfor

% for i in range(bnvars):
    ur[${i + bnvars}] = ul[${i + bnvars}];
% endfor
</%pyfr:macro>
