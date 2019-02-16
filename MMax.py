#!/usr/bin/python

import math
import numpy

from entangle import slater
from hv import energy_v
from htb import energy_tb

nr = 20

dr = float(2.0/(nr-1))

# r = chi/jt
chi = float(0.0)
r = [float(0.0) for i in range(nr)]

for i in range(nr):
    r[i] = float(i*dr)

jt = float(2.5)

gap = float(jt)
#gap = float(4*jt)
#gap = float(7*jt)

es0 = 0
et = -0.5*(gap - es0 - 4*jt)

n = 101

#calc. S for nearS_pts states nearly resonant with singlet
nearS_pts = 10

qvec = [int(0) for i in range(1,n-1,2)]
for i in range(1,n-1,2):
    qvec[(i-1)/2] = i
thq = numpy.multiply(float(2.0*math.pi/n),qvec)

#calculate entanglement

slat = [int(0) for j in range(nr*nearS_pts)]
slat = numpy.reshape(slat,(nr,nearS_pts))

e_tb = energy_tb(et,jt,thq,n)

h_tb = numpy.eye(len(thq))
for j in range(len(thq)):
    h_tb[j][j] = e_tb[j]
 
for i in range(nr):

    chi = jt*r[i]
    e_v = energy_v(chi,thq,n)
    h_v = numpy.reshape(e_v,(len(thq),len(thq)))

    energies, u = numpy.linalg.eig(h_tb + h_v)

    ut = numpy.conjugate(numpy.transpose(u))

    ind = 0
    etmp = 1.0e6
    for j in range(len(thq)):
        if etmp>abs(energies[j]):
            ind = j
            etmp = abs(energies[j])

    # points nearly resonant with singlet
    inds = [int(0) for j in range(nearS_pts)]
    for j in range(nearS_pts):
        inds[j] = ind-(nearS_pts/2)+j

    slat[i][:] = slater(u[:][inds],n,thq)

y = slat
ymax = int((n-1)/2)
y = y - 1
ymax = ymax - 1

for i in range(nr):
    # average over all nearS_pts calculated
    ymean[i] = float(numpy.mean(y[i][:])/ymax)

for i in range(nr):
    print("%f\t%f" % (r[i],ymean[i]))

