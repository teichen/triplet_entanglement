from math import pi
import numpy as np

from entangle import slater
from hv import energy_v
from htb import energy_tb

def MMax(r, n):
    """ calculate a quantum entanglement between two triplet excitations born out of
        singlet fission through a Slater rank for two-triplet points nearly resonant with
        the initial singlet excitation

        the governing two-triplet interaction strength, chi, is constrained by r = chi / jt,
        where chi is the triplet-triplet biexciton interaction

    Args:
        r (double): two-triplet separation
        n (int): size of the lattice

    Returns: ymean (np.array) mean Slater rank for triplet pairs nearly resonant 
             with initial singlet
    """
    jt = 2.5 # triplet-triplet resonance integral
    gap = jt # energy gap between two-triplet state and the ground state singlet

    es0 = 0 # the ground state singlet sets the zero
    et = -0.5 * (gap - es0 - 4 * jt)

    # calculate a Slater rank, S, for proximal states nearly resonant with singlet
    proximal_points = 10

    qvec = np.zeros(int((n-1) / 2)) # relative two-triplet momentum

    for idq in range(1, n-1, 2):
        qvec[int((idq-1) / 2)] = idq

    thq = 2.0 * pi / n * qvec

    # calculate entanglement
    nr = len(r)
    slat = np.zeros((nr, proximal_points))

    # calculate the tight-binding energy for increasing relative two-triplet momentum
    e_tb = energy_tb(et, jt, thq, n)

    h_tb = np.eye(len(thq))
    for idq in range(len(thq)):
        h_tb[idq][idq] = e_tb[idq]
     
    for idx in range(nr):

        chi = jt * r[idx]
        e_v = energy_v(chi, thq, n)
        h_v = np.reshape(e_v, (len(thq), len(thq)))

        # calculate the energies and unitary transformation
        energies, u = np.linalg.eig(h_tb + h_v)

        ut = np.conjugate(np.transpose(u))

        idx_nearest   = 0 # index nearest the singlet
        e_upper_bound = 1.0e6 # energy which bounds the singlet energy

        for idq in range(len(thq)):
            if e_upper_bound > abs(energies[idq]):
                idx_nearest = idq
                e_upper_bound = abs(energies[idq])

        # points nearly resonant with singlet
        idx_proximal = np.zeros(proximal_points, dtype=int)
        
        for idq in range(proximal_points):
            idx_proximal[idq] = idx_nearest - int(proximal_points / 2) + idq

        slat[idx, :] = slater(u[:, idx_proximal],n,thq)

    y    = slat
    ymax = int((n-1) / 2)
    y    = y - 1
    ymax = ymax - 1

    ymean = np.zeros(nr)

    for idx in range(nr):
        # average over all proximal_points calculated
        ymean[idx] = float(np.mean(y[idx, :]) / ymax)

    return ymean
