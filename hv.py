#!/usr/bin/python

from math import sin

def energy_v(chi, thq, n):
    """ biexciton interaction in the tight binding representation

    Args:
        chi (double): biexciton interaction energy between two triplet excitons
        thq (double): two triplet relative lattice momentum
        n (int): size of the lattice

    Returns: e_v (list) biexciton energy for a range of two triplet relative
             lattice momentum
    """

    e_v = [0.0 for j in range(len(thq)*len(thq))]

    for j in range(len(thq)):
        for k in range(len(thq)):
            e_v[j * len(thq) + k] = -(4 * chi / n) * sin(thq[j]) * sin(thq[k])

    return e_v
