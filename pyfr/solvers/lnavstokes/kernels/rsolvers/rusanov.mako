<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='lin_prim_normal_flux' params='du, ub, n, fn'>
    // ub = [rho_bar, u_bar..., p_bar]
    // du = [rho', u'..., p']
    fpdtype_t rhob = ub[0];
    fpdtype_t pb = ub[${nvars - 1}];

    fpdtype_t unb = ${pyfr.dot('n[{i}]', f'ub[{1} + {{i}}]', i=ndims)};
    fpdtype_t unp = ${pyfr.dot('n[{i}]', f'du[{1} + {{i}}]', i=ndims)};

    // rho' equation flux
    fn[0] = unb*du[0] + rhob*unp;

    // momentum/velocity equations flux
% for i in range(ndims):
    fn[${1 + i}] = unb*du[${1 + i}] + n[${i}]*du[${nvars - 1}]/rhob;
% endfor

    // pressure equation flux
    fn[${nvars - 1}] = unb*du[${nvars - 1}] + ${c['gamma']}*pb*unp;
</%pyfr:macro>


<%pyfr:macro name='rsolve' params='ul, ur, n, nf' externs='ub'>
    // ul/ur: primitive perturbations [rho', u', v'(,w'), p']
    // ub: primitive base states [rho, u, v(,w), p]

    fpdtype_t fnl[${nvars}], fnr[${nvars}];
    ${pyfr.expand('lin_prim_normal_flux', 'ul', 'ub', 'n', 'fnl')};
    ${pyfr.expand('lin_prim_normal_flux', 'ur', 'ub', 'n', 'fnr')};

    // Frozen Rusanov speed from base primitives
    fpdtype_t unb = ${pyfr.dot('n[{i}]', f'ub[{1} + {{i}}]', i=ndims)};
    fpdtype_t cl = sqrt(${c['gamma']}*ub[${nvars - 1}]/ub[0]);

    fpdtype_t a = 0.5*(fabs(unb) + cl);

% for i in range(nvars):
    nf[${i}] = 0.5*(fnl[${i}] + fnr[${i}]) + a*(ul[${i}] - ur[${i}]);
% endfor
</%pyfr:macro>
