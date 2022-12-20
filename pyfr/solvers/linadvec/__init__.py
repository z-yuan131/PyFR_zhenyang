# -*- coding: utf-8 -*-

from pyfr.solvers.linadvec.system import LinearAdvectionSystem
from pyfr.solvers.linadvec.elements import LinearAdvectionElements
from pyfr.solvers.linadvec.inters import (LinearAdvectionBCInters,
                                               LinearAdvectionIntInters,
                                               LinearAdvectionMPIInters)


"""
To calculate _basegrad_upts, nvars variables has been used from scal_upts which is too many.
Suggested to change it with slice or slicemat to reduce input to only Baseflow
variables. In memory allocation, use bvalloc rather than valloc. Some exploration
of the code is necessary.

Since in the interface, I used advecdiff inter, so shock-capturing somehow is needed
in the configuration file. This bug is needed to be eliminated.
"""
