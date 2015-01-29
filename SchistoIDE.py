__author__ = 'mario'

from numpy.fft import fft, ifft, fftshift
from numpy import array
import math
import matplotlib.pyplot as plt

# Variables definition
np = 256
x1 = 10.0
x0 = 0
dx = 2*x1/np
D = .01

mu = 0.0002
delta = 0.05
h = 100
f = 0.00005
b = .0005
s = 10
v = 0.0

#Note that the Critical Value is derived from the SchistosomiasisF.m in Matlab
Hcrit = 4.1368

c = h*f
nsteps = 2000
keyplot = [0, nsteps/2, nsteps-1]

#Discretization of initial distribution
x = []
for i in range(0,np-1):
    newrange = -x1 + i*dx
    x.append(newrange)
x = array(x)

I0 = array([0] * 255)

H0 = []
H0_element = array([abs(i) for i in x])
for i in H0_element:
    if i <=.25:
        H0_i = 20 * 1
        H0.append(H0_i)
    else:
        H0.append(0)
H0 = array(H0)

K = []

for i in x:
    equation = math.exp(-1*(i-v)**2/(4*D))/math.sqrt(4*math.pi*D)
    K.append(equation)

K_new = [i/(dx*sum(K)) for i in K]

fK = array(fft(K_new).T*dx)

#Graphing function starts here

critvalue = []
critvalueX = []
critvalueH = []
CritX = array(critvalueX)
CritH = array(critvalueH)


for j in range(0, nsteps):

    fI0 = array(fft(I0))
    I0 = array(fftshift(ifft(fK*fI0)).real)

    H = (1-mu)*H0+c*I0
    I = (1-delta)*I0+b*((s-I0)*H0**2)/(1+H0)

    for i in range(128, len(x)):
        if H[i] >= Hcrit:
            critvalue.append((i,j))
            critvalueH.append(j)
            critvalueX.append(i)

    if j in keyplot:
        plt.plot(x,H,'ro', linestyle = '-')
        plt.axis([x0,x1,-1,30])
        plt.xlabel('Distance')
        plt.ylabel('Worms (per host)')
        plt.title('Schistosomiasis: diffusion + advection')
        plt.show()

    I0 = I
    H0 = H

# print(critvalue)
# print(critvalueX)
# print(critvalueH)

Instantlistx = []
InstantlistH = []
for i in range(len(critvalueH)):
    if critvalueX[i] not in Instantlistx:
        Instantlistx.append(critvalueX[i])
        InstantlistH.append(critvalueH[i])

print((Instantlistx,InstantlistH))