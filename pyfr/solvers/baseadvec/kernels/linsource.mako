<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='linsource' params='u, linsrc' externs='ub, gb'>
    // Source terms for the linear solver
    fpdtype_t rhob = ub[0];
    fpdtype_t invrhob2 = 1.0/(rhob*rhob);

% if ndims == 2:
    fpdtype_t rhop = u[0], up = u[1], vp = u[2], pp = u[3];

    // Base gradients
    fpdtype_t ub_x = gb[0][1], ub_y = gb[1][1];
    fpdtype_t vb_x = gb[0][2], vb_y = gb[1][2];
    fpdtype_t pb_x = gb[0][3], pb_y = gb[1][3];

    // Contractions
    fpdtype_t up_grad_u1b = up*ub_x + vp*ub_y;     
    fpdtype_t up_grad_u2b = up*vb_x + vp*vb_y;     
    fpdtype_t div_ub = ub_x + vb_y;                
    fpdtype_t up_grad_pb = up*pb_x + vp*pb_y;      

    linsrc[1] += -up_grad_u1b + rhop*invrhob2*pb_x;
    linsrc[2] += -up_grad_u2b + rhop*invrhob2*pb_y;
    linsrc[3] += -(${c['gamma']} - 1.0)*(pp*div_ub + up_grad_pb);
% elif ndims == 3:
    fpdtype_t rhop = u[0], up = u[1], vp = u[2], wp = u[3], pp = u[4];

    // Base gradients
    fpdtype_t ub_x = gb[0][1], ub_y = gb[1][1], ub_z = gb[2][1];
    fpdtype_t vb_x = gb[0][2], vb_y = gb[1][2], vb_z = gb[2][2];
    fpdtype_t wb_x = gb[0][3], wb_y = gb[1][3], wb_z = gb[2][3];
    fpdtype_t pb_x = gb[0][4], pb_y = gb[1][4], pb_z = gb[2][4];

    // Contractions
    fpdtype_t up_grad_u1b = up*ub_x + vp*ub_y + wp*ub_z;   
    fpdtype_t up_grad_u2b = up*vb_x + vp*vb_y + wp*vb_z;   
    fpdtype_t up_grad_u3b = up*wb_x + vp*wb_y + wp*wb_z;   
    fpdtype_t div_ub = ub_x + vb_y + wb_z;                
    fpdtype_t up_grad_pb = up*pb_x + vp*pb_y + wp*pb_z;    

    linsrc[1] += -up_grad_u1b + rhop*invrhob2*pb_x;
    linsrc[2] += -up_grad_u2b + rhop*invrhob2*pb_y;
    linsrc[3] += -up_grad_u3b + rhop*invrhob2*pb_z;
    linsrc[4] += -(${c['gamma']} - 1.0)*(pp*div_ub + up_grad_pb);
% endif
</%pyfr:macro>
