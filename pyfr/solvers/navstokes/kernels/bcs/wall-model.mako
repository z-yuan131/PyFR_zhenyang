<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.common'/>

<%pyfr:macro name='bc_rsolve_state' params='ul, nl, ur' externs='ploc, t'>
    // Inviscid wall state for Riemann solve:
    // Reflect ONLY the normal component of velocity (slip wall).
    // u_r = u - 2*(u·n)*n

    fpdtype_t rho = ul[0];
    fpdtype_t invrho = 1.0/rho;

% if ndims == 2:
    fpdtype_t u0 = ul[1]*invrho;
    fpdtype_t u1 = ul[2]*invrho;

    fpdtype_t un = u0*nl[0] + u1*nl[1];

    fpdtype_t ur0 = u0 - 2.0*un*nl[0];
    fpdtype_t ur1 = u1 - 2.0*un*nl[1];

    ur[0] = rho;
    ur[1] = rho*ur0;
    ur[2] = rho*ur1;

    // Reflection preserves |u|, so keep total energy unchanged
    ur[3] = ul[3];

% elif ndims == 3:
    fpdtype_t u0 = ul[1]*invrho;
    fpdtype_t u1 = ul[2]*invrho;
    fpdtype_t u2 = ul[3]*invrho;

    fpdtype_t un = u0*nl[0] + u1*nl[1] + u2*nl[2];

    fpdtype_t ur0 = u0 - 2.0*un*nl[0];
    fpdtype_t ur1 = u1 - 2.0*un*nl[1];
    fpdtype_t ur2 = u2 - 2.0*un*nl[2];

    ur[0] = rho;
    ur[1] = rho*ur0;
    ur[2] = rho*ur1;
    ur[3] = rho*ur2;

    ur[4] = ul[4];
% endif
</%pyfr:macro>


<%pyfr:macro name='bc_ldg_state' params='ul, nl, ur' externs='ploc, t'>
    // WMLES boundary state for LDG/viscous gradients:
    // Enforce impermeability (u·n = 0) but allow tangential slip.
    // Work in conservative variables: ul = [rho, rhou, rhov, (rhow), E]

    fpdtype_t rho = ul[0];
    fpdtype_t invrho = 1.0/rho;

    // Velocity
% if ndims == 2:
    fpdtype_t u0 = ul[1]*invrho;
    fpdtype_t u1 = ul[2]*invrho;

    // Normal component
    fpdtype_t un = u0*nl[0] + u1*nl[1];

    // Remove normal component: u' = u - un*n
    fpdtype_t up0 = u0 - un*nl[0];
    fpdtype_t up1 = u1 - un*nl[1];

    // Conservative state
    ur[0] = rho;
    ur[1] = rho*up0;
    ur[2] = rho*up1;

    // Keep internal energy constant: E' = E - 0.5*rho*(un^2)
    ur[3] = ul[3] - 0.5*rho*un*un;

% elif ndims == 3:
    fpdtype_t u0 = ul[1]*invrho;
    fpdtype_t u1 = ul[2]*invrho;
    fpdtype_t u2 = ul[3]*invrho;

    fpdtype_t un = u0*nl[0] + u1*nl[1] + u2*nl[2];

    fpdtype_t up0 = u0 - un*nl[0];
    fpdtype_t up1 = u1 - un*nl[1];
    fpdtype_t up2 = u2 - un*nl[2];

    ur[0] = rho;
    ur[1] = rho*up0;
    ur[2] = rho*up1;
    ur[3] = rho*up2;

    ur[4] = ul[4] - 0.5*rho*un*un;
% endif
</%pyfr:macro>
