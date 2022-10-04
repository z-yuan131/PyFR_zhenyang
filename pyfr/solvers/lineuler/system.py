# -*- coding: utf-8 -*-

from pyfr.solvers.linadvec import LinearAdvectionSystem
from pyfr.solvers.lineuler.elements import LinearEulerElements
from pyfr.solvers.lineuler.inters import (LinearEulerIntInters, LinearEulerMPIInters,
                                       LinearEulerBaseBCInters)


class LinearEulerSystem(LinearAdvectionSystem):
    name = 'lineuler'

    elementscls = LinearEulerElements
    intinterscls = LinearEulerIntInters
    mpiinterscls = LinearEulerMPIInters
    bbcinterscls = LinearEulerBaseBCInters
