__author__ = 'mario'

from numpy.fft import fft, ifft, fftshift
from numpy import array, arange, meshgrid
import numpy
from scipy.stats import linregress
import math
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Variables definition
np = 256
x1 = 10.0
x0 = 0
dx = 2*x1/np
D = .01

mu = 0.0002
deltas = [0.01, .05, .1]
h = 100
f = 0.00005
b = .0005
s = 10
v = 0.0

#Note that the Critical Value is derived from the SchistosomiasisF.m in Matlab
Hcritlist = [.0776, .5645, 4.1368]
# Hcrit = 2.3

# Points at which the graph will draw
c = h*f
nsteps = 10000
keyplot = [0, nsteps-1]

#Graphing function starts here

critvalue = []
critvalueX = []
critvalueH = []
CritX = array(critvalueX)
CritH = array(critvalueH)

xlist = []
Hinitial = []
H1 = []
H2 = []

for delta in deltas:
    count = 0
    if delta == .01:
        Hcrit = Hcritlist[0]
    if delta == .05:
        Hcrit = Hcritlist[1]
        count = count + 1
    if delta == .1:
        Hcrit = Hcritlist[2]
        count = count + 1
    #Discretization of initial distribution
    x = []
    for i in range(0,np-1):
        newrange = -x1 + i*dx
        x.append(newrange)
        if count == 0:
            xlist.append(newrange)
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



    for j in range(0, nsteps):

        fI0 = array(fft(I0))
        I0 = array(fftshift(ifft(fK*fI0)).real)

        H = (1-mu)*H0+c*I0
        for i in H:
            if count == 0:
                Hinitial.append(i)
            if count == 1:
                H1.append(i)
            if count == 2:
                H2.append(i)
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
            # plt.show()

        I0 = I
        H0 = H
    Instantlistx = []
    InstantlistH = []
    for i in range(0, len(critvalueH)):
        if critvalueX[i] in Instantlistx:
            pass
        else:
            Instantlistx.append(critvalueX[i])
            InstantlistH.append(critvalueH[i])

    Instantlistx = [Instantlistx[i]-128 for i in range(len(Instantlistx))]
    # Instantlistx = array(Instantlistx)-128
    InstantlistH = (InstantlistH)
    print((Instantlistx,InstantlistH))

    plt.plot(Instantlistx, InstantlistH, 'b', linestyle = '-')
    plt.xlabel('Spatial Coordinate of Wave Front')
    plt.ylabel('Time to infection')
    plt.title('speed graph')
    plt.show()

    RegressX = Instantlistx[4:]
    RegressY = InstantlistH[4:]

    Regressionstats = linregress(RegressX,RegressY)
    print(Regressionstats)

#### 3D Plot ####
fig = plt.figure()
# ax = fig.gca(projection='3d')
ax = fig.add_subplot(111, projection='3d')

x = xlist
y = deltas

X,Y = meshgrid(x, y)
Z =array([Hinitial, H1, H2])

xs = [i for i in X]
ys = [i for i in Y]
zs = [i for i in Z]
c, m = 'r', 'o'
ax.scatter(xs, ys, zs, c=c, marker=m)
# ax.plot_surface(X, Z, rstride=8, cstride=8, alpha=0.3)
ax.plot_surface(X, Z, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(X, Y, Z, zdir='z', offset=-100, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
cset = ax.contour(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
print X,Y

def ThreeDPlot():

    x = y = arange(-3.0, 3.0, delta)
    X, Y = meshgrid(x, y)
    pass