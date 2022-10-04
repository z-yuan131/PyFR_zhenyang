# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='inviscid_flux' params='s, sb, f, p, v'>
    // Compute the average quantities
    fpdtype_t invrhob = 1.0/sb[0], Eb = sb[${nvars - 1}];

    fpdtype_t rhovb[${ndims}];
    fpdtype_t vb[${ndims}];
% for i in range(ndims):
    rhovb[${i}] = sb[${i + 1}];
    vb[${i}] = invrhob*sb[${i}]
% endfor

    // Averaged pressure
    pb = ${c['gamma'] - 1}*(Eb - 0.5*invrhob*${pyfr.dot('rhovb[{i}]', i=ndims)});





    fpdtype_t invrho = 1.0/s[0], E = s[${nvars - 1}];

    // Compute the velocities
    fpdtype_t rhov[${ndims}];
% for i in range(ndims):
    rhov[${i}] = s[${i + 1}];
    v[${i}] = invrho*rhov[${i}];
% endfor

    // Compute the pressure
    p = ${c['gamma'] - 1}*(E - 0.5*invrho*${pyfr.dot('rhov[{i}]', i=ndims)});


    // Density and energy fluxes
% for i in range(ndims):
    f[${i}][0] = sb[0]*v[${i}] + s[0]*vb[${i}];
    f[${i}][${nvars - 1}] = (Eb + pb)*v[${i}] + (E + p)*vb[${i}];
% endfor

    // Momentum fluxes
% for i, j in pyfr.ndrange(ndims, ndims):
    f[${i}][${j + 1}] = rhovb[${i}]*v[${j}]${' + p' if i == j else ''};
% endfor
</%pyfr:macro>
