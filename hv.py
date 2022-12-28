#!/usr/bin/python

import math

def energy_v(chi,thq,n):
    """ biexciton interaction in the tight binding representation
    """

    e_v = [0.0 for j in range(len(thq)*len(thq))]

    for j in range(len(thq)):
        for k in range(len(thq)):
            e_v[j * len(thq) + k] = -(4 * chi / n) * math.sin(thq[j]) * math.sin(thq[k])

    return e_v
