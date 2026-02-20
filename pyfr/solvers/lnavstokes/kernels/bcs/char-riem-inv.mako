<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.lnavstokes.kernels.bcs.common'/>

#include <stdio.h>


<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t, ub'>
    // Linear characteristic BC from solving the eigenvalue problem of a 
    // linear euler equation

    // printf("%f, %f, %f\n", ub[1], ul[1], ul[0]);


    fpdtype_t rhob = ub[0];
    fpdtype_t pb = ub[${nvars - 1}];
    fpdtype_t cb2 = ${c['gamma']}*pb/rhob;
    fpdtype_t cb = sqrt(cb2);
    fpdtype_t invrhobcb = 1.0/(rhob*cb);

    // Base normal velocity
    fpdtype_t unb = 0.0;
% for i in range(ndims):
    unb += nl[${i}]*ub[${i + 1}];
% endfor

    // Interior perturbation components
    fpdtype_t unr_l = 0.0;
% for i in range(ndims):
    unr_l += nl[${i}]*ul[${i + 1}];
% endfor
    fpdtype_t r_l = ul[0];
    fpdtype_t p_l = ul[${nvars - 1}];

    // Exterior/prescribed perturbation components
    fpdtype_t unr_e = 0.0;
% for i, v in enumerate('uvw'[:ndims]):
    fpdtype_t uext_${i} = ${c[v]};
    unr_e += nl[${i}]*uext_${i};
% endfor
    fpdtype_t r_e = ${c['rho']};
    fpdtype_t p_e = ${c['p']};

    // Acoustic invariants
    fpdtype_t wp_l = unr_l + p_l*invrhobcb;  // lambda = unb + cb
    fpdtype_t wm_l = unr_l - p_l*invrhobcb;  // lambda = unb - cb
    fpdtype_t wp_e = unr_e + p_e*invrhobcb;
    fpdtype_t wm_e = unr_e - p_e*invrhobcb;

    fpdtype_t wp_b = (unb + cb >= 0.0) ? wp_l : wp_e;
    fpdtype_t wm_b = (unb - cb >= 0.0) ? wm_l : wm_e;

    // Entropy-like convected invariant
    fpdtype_t s_l = p_l - cb2*r_l;           // lambda = unb
    fpdtype_t s_e = p_e - cb2*r_e;
    fpdtype_t s_b = (unb >= 0.0) ? s_l : s_e;

    // Reconstruct normal velocity, pressure, density perturbations
    fpdtype_t unr_b = 0.5*(wp_b + wm_b);
    fpdtype_t p_b = 0.5*rhob*cb*(wp_b - wm_b);
    fpdtype_t r_b = (p_b - s_b)/cb2;

    ur[0] = r_b;
    ur[${nvars - 1}] = p_b;

    // Tangential velocity perturbations are convected with speed unb
% for i in range(ndims):
    fpdtype_t ut_l_${i} = ul[${i + 1}] - unr_l*nl[${i}];
    fpdtype_t ut_e_${i} = uext_${i} - unr_e*nl[${i}];
    fpdtype_t ut_b_${i} = (unb >= 0.0) ? ut_l_${i} : ut_e_${i};
    ur[${i + 1}] = ut_b_${i} + unr_b*nl[${i}];
% endfor
</%pyfr:macro>

<%pyfr:alias name='bc_ldg_state' func='bc_rsolve_state'/>
<%pyfr:alias name='bc_ldg_grad_state' func='bc_common_grad_zero'/>
