__author__ = 'mario'

from numpy.fft import fft, ifft, fftshift
from numpy import array
from scipy.stats import linregress
import math
import matplotlib.pyplot as plt
from pylab import *

class Schistosomiasis():

    def __init__(self):
        # Variables definition
        self.np = np = 256
        self.x1 = x1 = 10.0
        self.x0 = 0
        self.dx = 2*x1/np
        self.D = .01

        self.mu = 0.0002
        self.delta = 0.05
        self.h = h = 100
        self.f = f = 0.00005
        self.b = .0005
        self.s = 10
        self.v = 0.05

        #Note that the Critical Value is derived from the SchistosomiasisF.m in Matlab
        self.Hcritlist = [.0776, .5645, 4.1368]
        # Hcrit = 2.3

        # Points at which the graph will draw
        self.chemical = h*f
        self.chemlist = [self.chemical, 2* self.chemical, 4 * self.chemical]
        self.nsteps = 10000
        self.keyplot = [0, self.nsteps-1]

        #Graphing function starts here

        self.critvalue = []
        self.critvalueX = []
        self.critvalueH = []
        self.CritX = array(self.critvalueX)
        self.CritH = array(self.critvalueH)

        self.MainLoop()

    def MainLoop(self):
        ### For Plotting Purposes, Lists not as important ###
        Count = 0
        xlist = []
        Hlist = []
        Hinitial = []
        H1 = []
        H2 = []
        H3 = []

        for c in self.chemlist:
            for delta in [0.01, .05, .1]:
                Count = Count + 1
                if delta == .01:
                    Hcrit = self.Hcritlist[0]
                if delta == .05:
                    Hcrit = self.Hcritlist[1]
                if delta == .1:
                    Hcrit = self.Hcritlist[2]
                #Discretization of initial distribution
                x = []
                for i in range(0,self.np-1):
                    newrange = -self.x1 + i*self.dx
                    x.append(newrange)
                x = array(x)
                xlist.append(x)

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
                    equation = math.exp(-1*(i-self.v)**2/(4*self.D))/math.sqrt(4*math.pi*self.D)
                    K.append(equation)

                K_new = [i/(self.dx*sum(K)) for i in K]

                fK = array(fft(K_new).T*self.dx)



                for j in range(0, self.nsteps):

                    fI0 = array(fft(I0))
                    I0 = array(fftshift(ifft(fK*fI0)).real)

                    H = (1-self.mu)*H0+c*I0
                    I = (1-delta)*I0+self.b*((self.s-I0)*H0**2)/(1+H0)

                    for i in range(128, len(x)):
                        if H[i] >= Hcrit:
                            self.critvalue.append((i,j))
                            self.critvalueH.append(j)
                            self.critvalueX.append(i)

                    if j in self.keyplot:
                        Hlist.append(H)
                        if Count == 1:
                            Hinitial.append(H)
                            print('Hi')
                        if Count == 2:
                            H1.append(H)
                            print('I am line 108, H1 append')
                        if Count == 3:
                            H2.append(H)
                            print('I am line 111')
                        if Count == 4:
                            H3.append(H)
                        else:
                            pass
                        # plt.plot(x,H,'ro', linestyle = '-')
                        # plt.axis([x0,x1,-1,30])
                        # plt.xlabel('Distance')
                        # plt.ylabel('Worms (per host)')
                        # plt.title('Schistosomiasis: diffusion + advection')
                        # plt.show()

                    I0 = I
                    H0 = H
                Instantlistx = []
                InstantlistH = []
                for i in range(0, len(self.critvalueH)):
                    if self.critvalueX[i] in Instantlistx:
                        pass
                    else:
                        Instantlistx.append(self.critvalueX[i])
                        InstantlistH.append(self.critvalueH[i])

                Instantlistx = [Instantlistx[i]-128 for i in range(len(Instantlistx))]
                # Instantlistx = array(Instantlistx)-128
                InstantlistH = (InstantlistH)
                print((Instantlistx,InstantlistH))



                plt.plot(Instantlistx, InstantlistH, 'b', linestyle = '-')
                plt.xlabel('Spatial Coordinate of Wave Front')
                plt.ylabel('Time to infection')
                plt.title('speed graph')
                if c == self.chemlist[0]:
                    if Count == 1:
                        savefig("Speed Graph Initial.png")
                    if Count == 2:
                        savefig("Speed Graph (.05).png")
                    if Count == 3:
                        savefig("Speed Graph (.1).png")
                if c == self.chemlist[1]:
                    if Count == 1:
                        savefig("Speed Graph Initial c2.png")
                    if Count == 2:
                        savefig("Speed Graph (.05) c2.png")
                    if Count == 3:
                        savefig("Speed Graph (.1) c2.png")
                if c == self.chemlist[2]:
                    if Count == 1:
                        savefig("Speed Graph Initial c3.png")
                    if Count == 2:
                        savefig("Speed Graph (.05) c3.png")
                    if Count == 3:
                        savefig("Speed Graph (.1) c3.png")
                # plt.show()

                RegressX = Instantlistx[4:]
                RegressY = InstantlistH[4:]

                Regressionstats = linregress(RegressX,RegressY)
                print(Regressionstats)



        # ### Plots of the distribution and spread of the Circuria ###
        # colors = ['r', 'g', 'b']
        # for i in range(len(xlist)):
        #     plt.figure()
        #     plt.plot(xlist[i], Hlist[i])
        #     plt.show()
        # #
        # plt.plot(xlist[0],Hinitial[0],'ro', linestyle = '-')
        # plt.plot(xlist[0],H1[0],'ro', linestyle = '-')
        # plt.plot(xlist[0],H2[0],'ro', linestyle = '-')
        # # plt.plot(x,H3[0],'ro', linestyle = '-')
        # plt.axis([x0,x1,-1,30])
        # plt.xlabel('Distance')
        # plt.ylabel('Worms (per host)')
        # plt.title('Schistosomiasis: diffusion + advection')
    # plt.show()

main = Schistosomiasis
if __name__ == '__main__':
    main()

