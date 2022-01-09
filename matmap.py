#!/usr/bin/python

import math
import numpy

def qr_householder(A):

    #returns V and tau such that H*A = alpha*e1
    #In this expression, e1 is the first column of eye(m), 
    #abs(alpha) = norm(A), and 
    #H = eye(m) - tau*V*V' is a Householder matrix.

    m = len(A)
    #n = len(A[0])
    Q = numpy.eye(m) 
    R = A.copy() 

    alpha = numpy.linalg.norm(A)
    alpha = -alpha*numpy.sign(A[0])

#    for j in range(n):
#        # Find H = I - tau*v*v' to put zeros below R[j,j]
#        x = R[j:, j]
#        normx = numpy.linalg.norm(x)
#        rho = -numpy.sign(x[0])
#        v1 = x[0] - rho * normx
#        v = x / v1
#        v[0] = 1
#        tau = -rho * v1 / normx

#        R[j:, :] = R[j:, :] - tau * numpy.outer(v, v).dot(R[j:, :])
#        Q[:, j:] = Q[:, j:] - tau * Q[:, j:].dot(numpy.outer(v, v))

#    return Q, R

    normR = numpy.linalg.norm(R)
    rho = -numpy.sign(R[0])
    v1 = R[0] - rho * normR
    v = R / v1
    v[0] = 1
    tau = -rho * v1 / normR

    return (v, tau, alpha)

def tridiag(A):

    #skew_tridiagonalize: Compute the tridiagonal form of A under unitary
    #congurence (orthogonal similarity, if A is real)

    #tridiag computes a skew-symmetric tridiagonal
    #matrix T and a unitary Q such that A = Q*T*transpose(Q)

    n = len(A)

    T = A

    Q = numpy.identity(n)

    for i in range(0,n-1):
    
        #Find a Householder vector to eliminate the i-th column

        Tvec = T[(i+1):n, i]
        (v, tau, alpha) = qr_householder(Tvec)

        T[i+1, i] = alpha
        T[i, i+1] = -alpha

        j = 2
        while i+j<n:
            T[i+j, i] = 0
            T[i, i+j] = 0
            j = j+1

        #Note: tau = 0 means the transformation is the identity

        if tau != 0.0:
            #Update the matrix block
            w = tau*numpy.matmul(T[(i+1):n, (i+1):n],numpy.conjugate(v))

            T[(i+1):n, (i+1):n] = T[(i+1):n, (i+1):n] + numpy.outer(v,w) - numpy.outer(w,v)

            #Accumulate the individual Householder reflections
            #Accumulate them in the form P_1*P_2*..., which is
            #(..*P_2*P_1)^dagger
            y = tau*numpy.matmul(Q[:, (i+1):n],v)
            Q[:, (i+1):n] = Q[:, (i+1):n] - numpy.outer(y,v)

    return (T, Q)

