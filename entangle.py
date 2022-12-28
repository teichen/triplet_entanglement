from math import sqrt, sin
import numpy

from matmap import tridiag

def slater(u, n, thq):
    """ slater rank
    """

    nq = len(u)
    proximal_points = len(u[0])

    slater_number = [int(0) for i in range(proximal_points)]

    for i in range(proximal_points):

        dmat = [0.0 for j in range(n*n)]
        dmat = numpy.reshape(dmat,(n,n))

        for j in range(1,n):
            for k in range(j+1,n+1):

                for l in range(nq):
                    dmat[j-1, k-1] = dmat[j-1, k-1] + u[l][i]*sin((k-j)*thq[l])

                dmat[j-1, k-1] = dmat[j-1, k-1]*(2/n)

                dmat[k-1, j-1] = -dmat[j-1, k-1]

        # now for indistinguishable fermions
        dmat = 0.5 * dmat / sqrt(numpy.trace(numpy.matmul(numpy.transpose(dmat), dmat))) # normalization

        # W. Wimmer
        # Algorithm 923: Efficient Numerical Computation of the Pfaffian 
        # for Dense and Banded Skew-Symmetric Matrices

        (dtilde, U) = tridiag(dmat)
        z = numpy.array(numpy.diagonal(dtilde, 1))

        for j in range(nq, n-1):
            z[j] = 0.0

        z2sum = 0.0
        for j in range(0, n-1):
            z2sum = z2sum + z[j] * z[j]

        for j in range(0, n-1):
            z[j] = z[j] * sqrt(0.25 / z2sum) # normalization    

        thrsh = float(1.0 / nq)

        stmp = 0
        for j in range(0,n-1):
            if z[j] > thrsh:
                stmp = stmp + 1

        slater_number[i] = stmp

        zoff = numpy.diagonal(dtilde, 2)

        z_null = 0
        for j in range(0, n-2):
            if zoff[j] == 0:
                z_null = z_null + 1
        if z_null != (n-2):
            slater_number[i] = 1e8

    return slater_number

