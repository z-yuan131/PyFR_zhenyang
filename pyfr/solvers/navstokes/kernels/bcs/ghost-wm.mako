<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.baseadvecdiff.kernels.artvisc'/>
<%include file='pyfr.solvers.euler.kernels.rsolvers.${rsolver}'/>
<%include file='pyfr.solvers.navstokes.kernels.flux'/>
<%include file='pyfr.solvers.navstokes.kernels.bcs.wmodel.${wmodel}'/>

<%pyfr:macro name='bc_common_flux_state' params='ul, gradul, artviscl, nl, magnl'>
    // Inviscid wall state (impermeable)
    fpdtype_t ur[${nvars}];
    ${pyfr.expand('bc_rsolve_state', 'ul', 'nl', 'ur')};

    // Inviscid Riemann solve
    fpdtype_t ficomm[${nvars}];
    ${pyfr.expand('rsolve', 'ul', 'ur', 'nl', 'ficomm')};

    // WMLES viscous flux: traction-only
    fpdtype_t fvr[${ndims}][${nvars}] = {{0}};

    // Compute wall shear stress vector tauw
    fpdtype_t tauw[${ndims}];
    ${pyfr.expand('compute_tau_wall', 'ul', 'gradul', 'nl', 'tauw')};

    // Fill viscous flux so that n · fvr gives tauw for momentum equations
% if ndims == 2:
    // rhou
    fvr[0][1] = tauw[0] * nl[0];
    fvr[1][1] = tauw[0] * nl[1];

    // rhov
    fvr[0][2] = tauw[1] * nl[0];
    fvr[1][2] = tauw[1] * nl[1];

    // Energy (adiabatic wall)
    fvr[0][3] = 0.0;
    fvr[1][3] = 0.0;
% elif ndims == 3:
    // rhou
    fvr[0][1] = tauw[0] * nl[0];
    fvr[1][1] = tauw[0] * nl[1];
    fvr[2][1] = tauw[0] * nl[2];

    // rhov
    fvr[0][2] = tauw[1] * nl[0];
    fvr[1][2] = tauw[1] * nl[1];
    fvr[2][2] = tauw[1] * nl[2];

    // rhow
    fvr[0][3] = tauw[2] * nl[0];
    fvr[1][3] = tauw[2] * nl[1];
    fvr[2][3] = tauw[2] * nl[2];

    // Energy (adiabatic)
    fvr[0][4] = 0.0;
    fvr[1][4] = 0.0;
    fvr[2][4] = 0.0;
% endif

    // Assemble common flux
    fpdtype_t fvcomm;

% for i in range(nvars):
    fvcomm = ${' + '.join(f'nl[{j}]*fvr[{j}][{i}]' for j in range(ndims))};
    ul[${i}] = magnl * (ficomm[${i}] + fvcomm);
% endfor
</%pyfr:macro>
