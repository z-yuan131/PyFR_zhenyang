<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

% if ndims == 2:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout' externs='ub, gb'>

    // Perturbation state (direct/LNS): [rho', u', v', p']
    fpdtype_t rho = uin[0], u = uin[1], v = uin[2], p = uin[3];

    // Base state: [rhob, ub, vb, pb]
    fpdtype_t rhob = ub[0], ubar = ub[1], vbar = ub[2], pb = ub[3];

    // Base reciprocals / constants
    fpdtype_t invrhob = 1.0/rhob;
    fpdtype_t invgmo  = 1.0/(${c['gamma']}-1.0);

    // Perturbation gradients
    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];

    fpdtype_t u_x = grad_uin[0][1];
    fpdtype_t u_y = grad_uin[1][1];
    fpdtype_t v_x = grad_uin[0][2];
    fpdtype_t v_y = grad_uin[1][2];

    fpdtype_t p_x = grad_uin[0][3];
    fpdtype_t p_y = grad_uin[1][3];

    // Base gradients 
    fpdtype_t rhob_x = gb[0][0];
    fpdtype_t rhob_y = gb[1][0];

    fpdtype_t ub_x = gb[0][1];
    fpdtype_t ub_y = gb[1][1];
    fpdtype_t vb_x = gb[0][2];
    fpdtype_t vb_y = gb[1][2];

    fpdtype_t pb_x = gb[0][3];
    fpdtype_t pb_y = gb[1][3];

% if visc_corr == 'sutherland':
    // Base cpT (= cp*T) using base p,rho: cpT = (gamma/(gamma-1)) * (pb/rhob)
    fpdtype_t cpT  = (${c['gamma']}/(${c['gamma']}-1.0))*(pb*invrhob);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;
    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}*Trat*sqrt(Trat)
                   / (cpT + ${c['cpTs']});

    // dmu/d(cpT) at base
    fpdtype_t dmudcpT = mu_c*(1.5/cpT - 1.0/(cpT + ${c['cpTs']}));

    // cpT' = cpT_b * (p'/pb - rho'/rhob)
    fpdtype_t mu_p = dmudcpT * cpT * (p/pb - rho*invrhob);
% else:
    fpdtype_t mu_c = ${c['mu']};
    fpdtype_t mu_p = 0.0;
% endif

    // Temperature variable here is Cv*T 
    // Tb = Cv*T_b = pb/(rhob*(gamma-1))
    fpdtype_t Tb   = invgmo*pb*invrhob;
    fpdtype_t Tb_x = invgmo*(invrhob*pb_x - pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y - pb*invrhob*invrhob*rhob_y);

    // T' (meaning Cv*T') = Tb*(p'/pb - rho'/rhob)
    // So grad(Cv*T') = Tb_x*(...) + Tb*grad(...)
    fpdtype_t T_x = Tb_x*(p/pb - rho*invrhob)
                    + Tb*( p_x/pb - p*pb_x/(pb*pb)
                    - rho_x*invrhob + rho*invrhob*invrhob*rhob_x );

    fpdtype_t T_y = Tb_y*(p/pb - rho*invrhob)
                    + Tb*( p_y/pb - p*pb_y/(pb*pb)
                    - rho_y*invrhob + rho*invrhob*invrhob*rhob_y );

    // Correct coefficients: (mu/rho)_b and (mu/rho)' 
    fpdtype_t muorho_b = mu_c*invrhob;

    // (mu/rho)' = mu'/rhob - mu_b*rho'/rhob^2
    fpdtype_t muorho_p = (mu_p - mu_c*(rho*invrhob))*invrhob;

    // Negated stresses used in viscous flux (consistent with your sign convention) 
    // From perturbation gradients using (mu/rho)_b
    fpdtype_t txx = -2.0*muorho_b*(u_x - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t tyy = -2.0*muorho_b*(v_y - ${1.0/3.0}*(u_x + v_y));
    fpdtype_t txy = -1.0*muorho_b*(v_x + u_y);

    // From base gradients using (mu/rho)'
    txx += -2.0*muorho_p*(ub_x - ${1.0/3.0}*(ub_x + vb_y));
    tyy += -2.0*muorho_p*(vb_y - ${1.0/3.0}*(ub_x + vb_y));
    txy += -1.0*muorho_p*(vb_x + ub_y);

    // Momentum viscous flux additions
    fout[0][1] += txx;
    fout[1][1] += txy;
    fout[0][2] += txy;
    fout[1][2] += tyy;

    // Base (negated) stresses for viscous work term in energy flux
    fpdtype_t tbxx = -2.0*muorho_b*(ub_x - ${1.0/3.0}*(ub_x + vb_y));
    fpdtype_t tbyy = -2.0*muorho_b*(vb_y - ${1.0/3.0}*(ub_x + vb_y));
    fpdtype_t tbxy = -1.0*muorho_b*(vb_x + ub_y);

    // Viscous work linearization: (tau u)' = tau'*ubar + tau_b*u'
    fout[0][3] += ubar*txx + vbar*txy + u*tbxx + v*tbxy;
    fout[1][3] += ubar*txy + vbar*tyy + u*tbxy + v*tbyy;

    // Heat conduction: -(k grad T)' = -k_b grad T' - k' grad Tb
    // with k = mu*cp/Pr and cp = gamma*cv => factor gamma/Pr since T here is cv*T
    fout[0][3] += -mu_c*${c['gamma']/c['Pr']}*T_x  - mu_p*${c['gamma']/c['Pr']}*Tb_x;
    fout[1][3] += -mu_c*${c['gamma']/c['Pr']}*T_y  - mu_p*${c['gamma']/c['Pr']}*Tb_y;

</%pyfr:macro>
% elif ndims == 3:
<%pyfr:macro name='viscous_flux_add' params='uin, grad_uin, fout' externs='ub, gb'>

    // Perturbation state: [rho', u', v', w', p']
    fpdtype_t rho = uin[0];
    fpdtype_t u   = uin[1];
    fpdtype_t v   = uin[2];
    fpdtype_t w   = uin[3];
    fpdtype_t p   = uin[4];

    // Base state
    fpdtype_t rhob = ub[0];
    fpdtype_t ubar = ub[1];
    fpdtype_t vbar = ub[2];
    fpdtype_t wbar = ub[3];
    fpdtype_t pb   = ub[4];

    fpdtype_t invrhob = 1.0/rhob;
    fpdtype_t invgmo  = 1.0/(${c['gamma']} - 1.0);

    // Perturbation gradients
    fpdtype_t rho_x = grad_uin[0][0];
    fpdtype_t rho_y = grad_uin[1][0];
    fpdtype_t rho_z = grad_uin[2][0];

    fpdtype_t u_x = grad_uin[0][1];
    fpdtype_t u_y = grad_uin[1][1];
    fpdtype_t u_z = grad_uin[2][1];

    fpdtype_t v_x = grad_uin[0][2];
    fpdtype_t v_y = grad_uin[1][2];
    fpdtype_t v_z = grad_uin[2][2];

    fpdtype_t w_x = grad_uin[0][3];
    fpdtype_t w_y = grad_uin[1][3];
    fpdtype_t w_z = grad_uin[2][3];

    fpdtype_t p_x = grad_uin[0][4];
    fpdtype_t p_y = grad_uin[1][4];
    fpdtype_t p_z = grad_uin[2][4];

    // Base gradients
    fpdtype_t rhob_x = gb[0][0];
    fpdtype_t rhob_y = gb[1][0];
    fpdtype_t rhob_z = gb[2][0];

    fpdtype_t ub_x = gb[0][1];
    fpdtype_t ub_y = gb[1][1];
    fpdtype_t ub_z = gb[2][1];

    fpdtype_t vb_x = gb[0][2];
    fpdtype_t vb_y = gb[1][2];
    fpdtype_t vb_z = gb[2][2];

    fpdtype_t wb_x = gb[0][3];
    fpdtype_t wb_y = gb[1][3];
    fpdtype_t wb_z = gb[2][3];

    fpdtype_t pb_x = gb[0][4];
    fpdtype_t pb_y = gb[1][4];
    fpdtype_t pb_z = gb[2][4];

% if visc_corr == 'sutherland':

    fpdtype_t cpT  = (${c['gamma']}/(${c['gamma']} - 1.0))*(pb*invrhob);
    fpdtype_t Trat = ${1/c['cpTref']}*cpT;

    fpdtype_t mu_c = ${c['mu']*(c['cpTref'] + c['cpTs'])}
                     * Trat*sqrt(Trat)
                     / (cpT + ${c['cpTs']});

    fpdtype_t dmudcpT = mu_c*(1.5/cpT - 1.0/(cpT + ${c['cpTs']}));
    fpdtype_t mu_p = dmudcpT * cpT * (p/pb - rho*invrhob);

% else:

    fpdtype_t mu_c = ${c['mu']};
    fpdtype_t mu_p = 0.0;

% endif

    // ---- Temperature (Cv*T) ----
    fpdtype_t Tb   = invgmo*pb*invrhob;

    fpdtype_t Tb_x = invgmo*(invrhob*pb_x - pb*invrhob*invrhob*rhob_x);
    fpdtype_t Tb_y = invgmo*(invrhob*pb_y - pb*invrhob*invrhob*rhob_y);
    fpdtype_t Tb_z = invgmo*(invrhob*pb_z - pb*invrhob*invrhob*rhob_z);

    fpdtype_t T_x = Tb_x*(p/pb - rho*invrhob)
                  + Tb*( p_x/pb - p*pb_x/(pb*pb)
                       - rho_x*invrhob + rho*invrhob*invrhob*rhob_x );

    fpdtype_t T_y = Tb_y*(p/pb - rho*invrhob)
                  + Tb*( p_y/pb - p*pb_y/(pb*pb)
                       - rho_y*invrhob + rho*invrhob*invrhob*rhob_y );

    fpdtype_t T_z = Tb_z*(p/pb - rho*invrhob)
                  + Tb*( p_z/pb - p*pb_z/(pb*pb)
                       - rho_z*invrhob + rho*invrhob*invrhob*rhob_z );

    // ---- (mu/rho) linearization ----
    fpdtype_t muorho_b = mu_c*invrhob;
    fpdtype_t muorho_p = (mu_p - mu_c*(rho*invrhob))*invrhob;

    fpdtype_t div_u  = u_x + v_y + w_z;
    fpdtype_t div_ub = ub_x + vb_y + wb_z;

    // ---- Stresses ----
    fpdtype_t txx = -2.0*muorho_b*(u_x - ${1.0/3.0}*div_u);
    fpdtype_t tyy = -2.0*muorho_b*(v_y - ${1.0/3.0}*div_u);
    fpdtype_t tzz = -2.0*muorho_b*(w_z - ${1.0/3.0}*div_u);

    fpdtype_t txy = -muorho_b*(v_x + u_y);
    fpdtype_t txz = -muorho_b*(w_x + u_z);
    fpdtype_t tyz = -muorho_b*(w_y + v_z);

    txx += -2.0*muorho_p*(ub_x - ${1.0/3.0}*div_ub);
    tyy += -2.0*muorho_p*(vb_y - ${1.0/3.0}*div_ub);
    tzz += -2.0*muorho_p*(wb_z - ${1.0/3.0}*div_ub);

    txy += -muorho_p*(vb_x + ub_y);
    txz += -muorho_p*(wb_x + ub_z);
    tyz += -muorho_p*(wb_y + vb_z);

    // Momentum flux
    fout[0][1] += txx;  fout[1][1] += txy;  fout[2][1] += txz;
    fout[0][2] += txy;  fout[1][2] += tyy;  fout[2][2] += tyz;
    fout[0][3] += txz;  fout[1][3] += tyz;  fout[2][3] += tzz;

    // ---- Base stresses for energy viscous work ----
    fpdtype_t tbxx = -2.0*muorho_b*(ub_x - ${1.0/3.0}*div_ub);
    fpdtype_t tbyy = -2.0*muorho_b*(vb_y - ${1.0/3.0}*div_ub);
    fpdtype_t tbzz = -2.0*muorho_b*(wb_z - ${1.0/3.0}*div_ub);

    fpdtype_t tbxy = -muorho_b*(vb_x + ub_y);
    fpdtype_t tbxz = -muorho_b*(wb_x + ub_z);
    fpdtype_t tbyz = -muorho_b*(wb_y + vb_z);

    // Energy viscous work
    fout[0][4] += ubar*txx + vbar*txy + wbar*txz
                + u*tbxx   + v*tbxy   + w*tbxz;

    fout[1][4] += ubar*txy + vbar*tyy + wbar*tyz
                + u*tbxy   + v*tbyy   + w*tbyz;

    fout[2][4] += ubar*txz + vbar*tyz + wbar*tzz
                + u*tbxz   + v*tbyz   + w*tbzz;

    // Heat conduction
    fout[0][4] += -mu_c*${c['gamma']/c['Pr']}*T_x  - mu_p*${c['gamma']/c['Pr']}*Tb_x;
    fout[1][4] += -mu_c*${c['gamma']/c['Pr']}*T_y  - mu_p*${c['gamma']/c['Pr']}*Tb_y;
    fout[2][4] += -mu_c*${c['gamma']/c['Pr']}*T_z  - mu_p*${c['gamma']/c['Pr']}*Tb_z;

</%pyfr:macro>
% endif
