# -*- coding: utf-8 -*-
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>
<%include file='pyfr.solvers.lineuler.kernels.flux'/>

<%pyfr:macro name='rsolve' params='ul, ur, ubl, ubr, n, nf'>
    // Compute the left and right fluxes + velocities and pressures
    fpdtype_t fl[${ndims}][${nvars}], fr[${ndims}][${nvars}];
    fpdtype_t vl[${ndims}], vr[${ndims}];
    fpdtype_t pl, pr;
    fpdtype_t pbl, pbr;


    ${pyfr.expand('inviscid_flux', 'ul', 'ubl', 'fl', 'pl', 'vl')};
    ${pyfr.expand('inviscid_flux', 'ur', 'ubr', 'fr', 'pr', 'vr')};


    // Sum the left and right velocities and take the normal
    fpdtype_t nv = ${pyfr.dot('n[{i}]', 'vl[{i}] + vr[{i}]', i=ndims)};

    // Get baseflow pressures
    pbl = ubl[${nvars - 1}];
    pbr = ubr[${nvars - 1}];

    // Estimate the maximum wave speed / 2
    fpdtype_t a = sqrt(${0.25*c['gamma']}*(pbl + pbr + pl + pr)/(ubl[0] + ubr[0] + ul[0] + ur[0]))
                + 0.25*fabs(nv);

    // Output
% for i in range(nvars):
    nf[${i}] = 0.5*(${' + '.join('n[{j}]*(fl[{j}][{i}] + fr[{j}][{i}])'
                                 .format(i=i, j=j) for j in range(ndims))})
             + a*(ul[${i}] - ur[${i}]);
% endfor




</%pyfr:macro>
