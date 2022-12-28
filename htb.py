#!/usr/bin/python

import math

def energy_tb(et, jt, thq, n):
    """ tight binding energy
    """

    e_tb = [0.0 for j in range(len(thq))]

    for j in range(len(thq)):
        e_tb[j] = 2 * et - 4 * jt * math.cos(thq[j])
    
    return e_tb

