<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

#include <stdio.h>

<%pyfr:macro name='compute_tau_wall' params='ul, gradul, nl, tauw' externs='nsoln, ym, tauw_exp'>
    // Constants
    fpdtype_t kappa = 0.41;
    fpdtype_t B = 5.2;
    fpdtype_t eps = 1e-12;

    // Extract primitive variables 
    fpdtype_t rho = ul[0];
    fpdtype_t rcprho = 1.0/rho;
    fpdtype_t u   = rcprho * ul[1];
    fpdtype_t v   = rcprho * ul[2];

% if ndims == 3:
    fpdtype_t w   = rcprho * ul[3];
    fpdtype_t E   = ul[4];
% else:
    fpdtype_t E   = ul[3];
% endif

    // Velocity vector
    fpdtype_t vel[${ndims}] = {
        u,
        v
% if ndims == 3:
        , w
% endif
    };

    // Molecular viscosity at wall
% if visc_corr == 'sutherland':
% if ndims == 3:
    fpdtype_t ke = 0.5*(uw*uw + vw*vw + ww*ww);
% else:
    fpdtype_t ke = 0.5*(uw*uw + vw*vw);
% endif
    // Compute the temperature and viscosity
    fpdtype_t cpT = ${c['gamma']}*(rcprho*E - ke);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}
                   * Trat * sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif
    fpdtype_t nu = mu_c / rho;



    // Approximate sampling state at the edge of wall model 
    
    fpdtype_t velm[${ndims}] = {
        nsoln[1], nsoln[2]
% if ndims == 3:
        , nsoln[3]
% endif
    };

    // Remove normal component
    fpdtype_t unm = 0.0;
% for d in range(ndims):
    unm += velm[${d}] * nl[${d}];
% endfor

    fpdtype_t utm[${ndims}];
    fpdtype_t Ut2 = 0.0;
% for d in range(ndims):
    utm[${d}] = velm[${d}] - unm * nl[${d}];
    Ut2 += utm[${d}] * utm[${d}];
% endfor
    fpdtype_t Ut = sqrt(Ut2);


    // Solve log-law for friction velocity
    fpdtype_t utau = 0.0;

    if (Ut > 1e-10)
    {
        utau = sqrt(nu * Ut / ym + eps);

        for (int it = 0; it < 6; ++it)
        {
            fpdtype_t yplus = ym * utau / nu;
            fpdtype_t f  = Ut / (utau + eps)
                         - (1.0/kappa) * log(yplus + eps)
                         - B;

            fpdtype_t df = -Ut / ((utau + eps)*(utau + eps))
                           - (1.0 / kappa) / (utau + eps);

            utau -= f / (df + eps);
            utau = fmax(utau, eps);
        }
    }


    // printf("%f", utau);
    // printf("%f\n", mu_c);

    // Wall shear stress magnitude
    fpdtype_t tau_mag = rho * utau * utau;

    // Wall shear stress
% for i in range(ndims):
    tauw[${i}] = tau_mag * utm[${i}] / (Ut + eps);
% endfor

    // Output the wall shear stress magnitude
    // stressinfo
    tauw_exp[0] = ym;
    tauw_exp[1] = utau;
    tauw_exp[2] = nu;

</%pyfr:macro>
