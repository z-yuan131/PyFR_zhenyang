<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

#include <stdio.h>

<%pyfr:macro name='compute_tau_wall' params='ul, gradul, nl, tauw' externs='nsoln, ym, tauw_exp'>
    // ---------------- Constants ----------------
    fpdtype_t kappa = 0.41;
    fpdtype_t eps   = 1e-12;

    // ---------------- Wall/interior-face state ("wall reference") ----------------
    fpdtype_t rho_w  = ul[0];
    fpdtype_t rcprho = 1.0/(rho_w + eps);

    fpdtype_t u = rcprho * ul[1];
    fpdtype_t v = rcprho * ul[2];
% if ndims == 3:
    fpdtype_t w = rcprho * ul[3];
    fpdtype_t E = ul[4];
% else:
    fpdtype_t E = ul[3];
% endif

    // Molecular viscosity at wall reference
% if visc_corr == 'sutherland':
% if ndims == 3:
    fpdtype_t ke = 0.5*(u*u + v*v + w*w);
% else:
    fpdtype_t ke = 0.5*(u*u + v*v);
% endif
    fpdtype_t cpT  = ${c['gamma']}*(rcprho*E - ke);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_w = ${c['mu']*(c['cpTref'] + c['cpTs'])}
                   * Trat * sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_w = ${c['mu']};
% endif

    fpdtype_t nu_w = mu_w / (rho_w + eps);

    // ---------------- Sampled state at y = ym ----------------
    fpdtype_t rho_m = nsoln[0];

    fpdtype_t velm[${ndims}] = {
        nsoln[1], nsoln[2]
% if ndims == 3:
        , nsoln[3]
% endif
    };

    // Remove normal component -> tangential velocity at ym
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

    // ---------------- van Driest transform (simple single-point form) ----------------
    // U_vd ≈ U * sqrt(rho_m/rho_w)
    fpdtype_t sqrt_rhor = sqrt((rho_m + eps)/(rho_w + eps));
    fpdtype_t Ut_vd = Ut * sqrt_rhor;

    // ---------------- Spalding law: solve for utau ----------------
    // Define:
    //   u+  = U_vd / utau
    //   y+  = ym * utau / nu_w
    // Spalding implicit relation:
    //   y+ = u+ + (1/kappa) [ exp(kappa u+) - 1 - kappa u+ - (kappa u+)^2/2 - (kappa u+)^3/6 ]
    //
    // Unknown: utau (>=0). Use Newton on F(utau)=0.

    fpdtype_t utau = 0.0;

    if (Ut_vd > 1e-10)
    {
        // Initial guess: log-law-ish
        utau = sqrt(nu_w * Ut_vd / (ym + eps));

        for (int it = 0; it < 8; ++it)
        {
            fpdtype_t uplus = Ut_vd / (utau + eps);
            fpdtype_t yplus = ym * utau / nu_w;

            // To avoid overflow in exp for very large uplus:
            // clamp kappa*uplus to a reasonable range
            fpdtype_t ku = kappa * uplus;
            ku = fmin(ku, (fpdtype_t)50.0);

            fpdtype_t expku = exp(ku);

            // Spalding RHS(u+) = u+ + 1/kappa * (exp(ku) - 1 - ku - ku^2/2 - ku^3/6)
            fpdtype_t ku2 = ku*ku;
            fpdtype_t ku3 = ku2*ku;

            fpdtype_t RHS = uplus + (1.0/kappa) * (expku - 1.0 - ku - 0.5*ku2 - (1.0/6.0)*ku3);

            // Residual F = RHS - y+
            fpdtype_t F = RHS - yplus;

            // dRHS/du+:
            // RHS = u+ + 1/kappa * (exp(ku) - 1 - ku - ku^2/2 - ku^3/6)
            // d/du+ [exp(ku)] = exp(ku) * kappa
            // d/du+ [-ku]     = -kappa
            // d/du+ [-ku^2/2] = -(ku)*kappa
            // d/du+ [-ku^3/6] = -(ku^2/2)*kappa
            //
            // => dRHS/du+ = 1 + 1/kappa * [kappa*exp(ku) - kappa - kappa*ku - kappa*ku^2/2]
            //            = 1 + (exp(ku) - 1 - ku - ku^2/2)
            fpdtype_t dRHS_du = 1.0 + (expku - 1.0 - ku - 0.5*ku2);

            // du+/dutau = d/dutau (U_vd/utau) = -U_vd/utau^2 = -(u+)/utau
            fpdtype_t duplus_dutau = -Ut_vd / ((utau + eps)*(utau + eps));

            // dy+/dutau = ym/nu_w
            fpdtype_t dyplus_dutau = ym / (nu_w + eps);

            // dF/dutau = dRHS/du+ * du+/dutau - dy+/dutau
            fpdtype_t dF = dRHS_du * duplus_dutau - dyplus_dutau;

            // Newton step
            fpdtype_t dut = -F / (dF + eps);

            // Damping for robustness (optional)
            // Limit relative update to avoid overshoot
            fpdtype_t maxrel = 0.5;
            fpdtype_t lim = maxrel * (utau + eps);
            dut = fmax(fmin(dut, lim), -lim);

            utau += dut;
            utau = fmax(utau, eps);

            // Early exit
            if (fabs(F) < 1e-10)
                break;
        }
    }

    // ---------------- Wall shear stress vector ----------------
    fpdtype_t tau_mag = (rho_w + eps) * utau * utau;

% for i in range(ndims):
    tauw[${i}] = tau_mag * utm[${i}] / (Ut + eps);
% endfor

    // Optional debug output
    tauw_exp[0] = ym;
    tauw_exp[1] = utau;
    tauw_exp[2] = nu_w;

</%pyfr:macro>
