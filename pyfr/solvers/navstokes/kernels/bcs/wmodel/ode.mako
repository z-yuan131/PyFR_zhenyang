<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='compute_tau_wall' params='ul, gradul, nl, tauw' externs='nsoln, ym, tauw_exp'>

    // ---------------- Constants ----------------
    fpdtype_t kappa = 0.41;
    fpdtype_t Aplus = 17.0;
    fpdtype_t eps   = 1e-12;

    fpdtype_t gamma = ${c['gamma']};

    enum { NWM = 16 };
    fpdtype_t dy = ym / (NWM - 1);

    // ---------------- Wall primitive state ----------------
    fpdtype_t rho    = ul[0];
    fpdtype_t rcprho = 1.0 / (rho + eps);

    fpdtype_t uw = rcprho * ul[1];
    fpdtype_t vw = rcprho * ul[2];
% if ndims == 3:
    fpdtype_t ww = rcprho * ul[3];
    fpdtype_t Ew = ul[4];
% else:
    fpdtype_t Ew = ul[3];
% endif

    // ---------------- Sample at y = ym ----------------
    fpdtype_t velm[${ndims}] = {
        nsoln[1], nsoln[2]
% if ndims == 3:
        , nsoln[3]
% endif
    };

    // Tangential velocity
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
    fpdtype_t Ut = sqrt(fmax(Ut2, 0.0));

    // ---------------- Molecular viscosity ----------------
% if visc_corr == 'sutherland':
% if ndims == 3:
    fpdtype_t ke = 0.5*(uw*uw + vw*vw + ww*ww);
% else:
    fpdtype_t ke = 0.5*(uw*uw + vw*vw);
% endif
    fpdtype_t cpT  = gamma*(rcprho*Ew - ke);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}
                   * Trat * sqrt(Trat)
                   / (cpT + ${c['cpTs']});
% else:
    fpdtype_t mu_c = ${c['mu']};
% endif

    fpdtype_t nuw = mu_c / rho;

    // ---------------- ODE wall model ----------------
    fpdtype_t utau = 0.0;

    if (Ut > 1e-10)
    {
        utau = sqrt(fmax(nuw * Ut / (ym + eps), eps));

        fpdtype_t y[NWM], u[NWM], nut[NWM];

        for (int i = 0; i < NWM; ++i)
        {
            y[i] = i * dy;
            u[i] = Ut * y[i] / (ym + eps);
        }

        for (int it = 0; it < 8; ++it)
        {
            for (int i = 0; i < NWM; ++i)
            {
                fpdtype_t yplus = y[i] * utau / (nuw + eps);
                fpdtype_t fd = 1.0 - exp(-yplus / Aplus);
                nut[i] = kappa * y[i] * utau * fd * fd;
            }

            u[0] = 0.0;
            u[NWM-1] = Ut;

            for (int sweep = 0; sweep < 8; ++sweep)
            {
                for (int i = 1; i < NWM-1; ++i)
                {
                    fpdtype_t nu_p = nuw + 0.5*(nut[i] + nut[i+1]);
                    fpdtype_t nu_m = nuw + 0.5*(nut[i] + nut[i-1]);
                    u[i] = (nu_p*u[i+1] + nu_m*u[i-1]) /
                           (nu_p + nu_m + eps);
                }
            }

            fpdtype_t dudy = (u[1] - u[0]) / (dy + eps);
            fpdtype_t tau_eff = fmax((nuw + nut[1]) * fmax(dudy, 0.0), eps);
            utau = 0.7*utau + 0.3*sqrt(tau_eff);
        }
    }

    // ---------------- Wall shear stress ----------------
    fpdtype_t tau_mag = rho * utau * utau;

% for d in range(ndims):
    tauw[${d}] = tau_mag * utm[${d}] / (Ut + eps);
% endfor
    

    // stressinfo
    tauw_exp[0] = ym;
    tauw_exp[1] = utau;
    tauw_exp[2] = nuw;

</%pyfr:macro>
