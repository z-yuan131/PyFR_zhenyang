# -*- coding: utf-8 -*-

from pyfr.solvers.linadvec import LinearAdvectionSystem
from pyfr.solvers.lineuler.elements import LinearEulerElements
from pyfr.solvers.lineuler.inters import (LinearEulerIntInters,
                    LinearEulerMPIInters, LinearEulerBaseBCInters)


class LinearEulerSystem(LinearAdvectionSystem):
    name = 'linear-euler'

    elementscls = LinearEulerElements
    intinterscls = LinearEulerIntInters
    mpiinterscls = LinearEulerMPIInters
    bbcinterscls = LinearEulerBaseBCInters
