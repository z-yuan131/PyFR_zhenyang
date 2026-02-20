<%inherit file='base'/>
<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%include file='pyfr.solvers.lnavstokes.kernels.bcs.${bctype}'/>

% if bccfluxstate:
<%include file='pyfr.solvers.lnavstokes.kernels.bcs.${bccfluxstate}'/>
% endif

#include <stdio.h>


<%pyfr:kernel name='bccflux' ndim='1'
              ul='inout view fpdtype_t[${str(nvars)}]'
              gradul='in view fpdtype_t[${str(ndims)}][${str(nvars)}]'
              artviscl='in view fpdtype_t'
              nl='in fpdtype_t[${str(ndims)}]'>
    fpdtype_t mag_nl = sqrt(${pyfr.dot('nl[{i}]', i=ndims)});
    fpdtype_t norm_nl[] = ${pyfr.array('(1 / mag_nl)*nl[{i}]', i=ndims)};


    // printf("%f, %f, %f, %f\n", gradul[0][0], gradul[0][3], gradul[0][2], gradul[0][1]);


    ${pyfr.expand('bc_common_flux_state', 'ul', 'gradul', 'artviscl', 'norm_nl', 'mag_nl')};
</%pyfr:kernel>
