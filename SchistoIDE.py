__author__ = 'mario'

from numpy.fft import fft, ifft, fftshift
from numpy import array
import math
import matplotlib.pyplot as plt

# Variables definition
np = 256
x1 = 10
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

c = h*f
nsteps = 10

#Discretization of initial distribution
x = []
for i in range(0,np-1):
    newrange = -x1 + i*dx
    x.append(newrange)

I0 = array([0] * 255)

H0 = array([abs(i) for i in x])
print(H0)
# H0
K = []

for i in x:
    equation = math.exp(-1*(i-v)**2/(4*D))/math.sqrt(4*math.pi*D)
    K.append(equation)

K_new = [i/(dx*sum(K)) for i in K]

fK = array(fft(K_new).T*dx)

#Graphing function starts here

for j in range(0, nsteps):

    fI0 = array(fft(I0))
    I0 = array(fftshift(ifft(fK*fI0)).real)

    H = (1-mu)*H0+c*I0
    I = (1-delta)*I0+b*((s-I0)*H0**2)/(1+H0)

    plt.plot(x,H,'ro')
    plt.axis([0,10,0,100])
    plt.show()
    plt.xlabel('distance')
    plt.ylabel('Probability')
    plt.title('Schistosomiasis: diffusion + advection')

    I0_new = I
    H0 = H


    pass

    # for i in range(len(fK)):


