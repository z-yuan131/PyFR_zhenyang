<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.lnavstokes.kernels.bcs.common'/>

#include <stdio.h>


<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = -ul[${i + 1}];
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
</%pyfr:macro>

<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t, ub'>
    ur[0] = ul[0];
% for i in range(ndims):
    ur[${i + 1}] = 0.0;
% endfor
    ur[${nvars - 1}] = ul[${nvars - 1}];
    // If Iso-thermal ur[${nvars - 1}] = ubl[${nvars - 1}]*ul[0]/ubl[0];


    // printf("%f, %f, %f, %f, %f\n", ub[1], ul[0], ul[1], ul[2], ul[3]);

</%pyfr:macro>



<%pyfr:macro name='bc_ldg_grad_state' params='ur, nl, grad_ul, grad_ur' externs='ub, gb'>
    // Copy all fluid-side gradients across to wall-side gradients
    ${pyfr.expand('bc_common_grad_copy', 'ur', 'nl', 'grad_ul', 'grad_ur')};

    // Enforce adiabatic wall in linear primitive form: n . grad(Cv*T') = 0
    fpdtype_t rho = ur[0], p = ur[${nvars - 1}];
    fpdtype_t rhob = ub[0], pb = ub[${nvars - 1}];
    fpdtype_t invrhob = 1.0/rhob;
    fpdtype_t invgmo = 1.0/(${c['gamma']} - 1.0);
    fpdtype_t Tb = invgmo*pb*invrhob;
    fpdtype_t phi = p/pb - rho*invrhob;
    fpdtype_t nTl = 0.0;

% if ndims == 2:
    fpdtype_t rhob_x = gb[0][0], rhob_y = gb[1][0];
    fpdtype_t pb_x = gb[0][${nvars - 1}], pb_y = gb[1][${nvars - 1}];
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x - pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y - pb*invrhob*invrhob*rhob_y);

    fpdtype_t rho_x = grad_ul[0][0], rho_y = grad_ul[1][0];
    fpdtype_t p_x = grad_ul[0][${nvars - 1}], p_y = grad_ul[1][${nvars - 1}];

    fpdtype_t Tl_x = Tb_x*phi
                   + Tb*(p_x/pb - p*pb_x/(pb*pb)
                   - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t Tl_y = Tb_y*phi
                   + Tb*(p_y/pb - p*pb_y/(pb*pb)
                   - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);

    nTl = nl[0]*Tl_x + nl[1]*Tl_y;

    // fpdtype_t bb = Tb*(p_x/pb - p*pb_x/(pb*pb));

    // printf("%f, %f, %f, %f\n", nTl, grad_ul[0][3], grad_ul[0][2], grad_ul[0][1]);

% elif ndims == 3:
    fpdtype_t rhob_x = gb[0][0], rhob_y = gb[1][0], rhob_z = gb[2][0];
    fpdtype_t pb_x = gb[0][${nvars - 1}], pb_y = gb[1][${nvars - 1}], pb_z = gb[2][${nvars - 1}];
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x - pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y - pb*invrhob*invrhob*rhob_y);
    fpdtype_t Tb_z = invgmo*(invrhob*pb_z - pb*invrhob*invrhob*rhob_z);

    fpdtype_t rho_x = grad_ul[0][0], rho_y = grad_ul[1][0], rho_z = grad_ul[2][0];
    fpdtype_t p_x = grad_ul[0][${nvars - 1}], p_y = grad_ul[1][${nvars - 1}], p_z = grad_ul[2][${nvars - 1}];

    fpdtype_t Tl_x = Tb_x*phi
                   + Tb*(p_x/pb - p*pb_x/(pb*pb)
                   - rho_x*invrhob + rho*invrhob*invrhob*rhob_x);
    fpdtype_t Tl_y = Tb_y*phi
                   + Tb*(p_y/pb - p*pb_y/(pb*pb)
                   - rho_y*invrhob + rho*invrhob*invrhob*rhob_y);
    fpdtype_t Tl_z = Tb_z*phi
                   + Tb*(p_z/pb - p*pb_z/(pb*pb)
                   - rho_z*invrhob + rho*invrhob*invrhob*rhob_z);

    nTl = nl[0]*Tl_x + nl[1]*Tl_y + nl[2]*Tl_z;
% endif

    // T' depends linearly on p' gradient with coefficient Tb/pb.
    // Correct only grad(p') to remove the normal component of grad(T').
    fpdtype_t coef = pb/Tb;
% for i in range(ndims):
    grad_ur[${i}][${nvars - 1}] -= coef*nl[${i}]*nTl;
% endfor
</%pyfr:macro>
