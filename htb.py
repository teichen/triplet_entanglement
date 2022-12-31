#!/usr/bin/python

from math import cos

def energy_tb(et, jt, thq, n):
    """ tight binding energy

    Args:
        et (double): single triplet local electronic energy
        jt (double): triplet resonance integral
        thq (double): two triplet relative lattice momentum

    Returns: e_tb (list) tight binding energies for a range of
             two triplet relative lattice momentum
    """

    e_tb = [0.0 for j in range(len(thq))]

    for j in range(len(thq)):
        e_tb[j] = 2 * et - 4 * jt * cos(thq[j])
    
    return e_tb

