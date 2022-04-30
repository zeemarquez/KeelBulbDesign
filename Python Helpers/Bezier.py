from msilib.schema import Class
import numpy as np

class BezierCurve:
    def __init__(self,c0,z0,c1,z1,z2,c3,z3,I):
        self.c0 = c0
        self.z0 = z0
        self.c1 = c1
        self.z1 = z1
        self.z2 = z2
        self.c3 = c3
        self.z3 = z3
        self.I = I

        self.c2 = self.get_c2()
    
    def Bezier(self,t):

        c = ((1-t)**3)*self.c0 + 3*((1-t)**2)*t*self.c1 + 3*((1-t)*t**2)*self.c2 + self.c3*t**3
        z = ((1-t)**3)*self.z0 + 3*((1-t)**2)*t*self.z1 + 3*((1-t)*t**2)*self.z2 + self.z3*t**3
        return (c, z)

    def getBezierCurve(self, resolution = 0.0001):

        c_tab , z_tab = [self.c0],[self.z0]

        for t in np.arange(0,1,resolution):
            c,z = self.Bezier(t)
            c_tab.append(c)
            z_tab.append(z)

        c_tab.append(self.c3)
        z_tab.append(self.z3)

        self.bezierCurve = (c_tab,z_tab)


        return (c_tab,z_tab)

    def get_c2(self):
        c0 = self.c0
        z0 = self.z0
        c1 = self.c1
        z1 = self.z1
        z2 = self.z2
        c3 = self.c3
        z3 = self.z3
        I = self.I

        sqrt = lambda x: x**0.5
        solution = (-10*c0*z0 + 6*c0*z2 + 4*c0*z3 - 15*c1*z0 - 9*c1*z1 + 9*c1*z2 + 15*c1*z3 - 5*c3*z0 - 15*c3*z1 - 15*c3*z2 + 35*c3*z3 - sqrt(3)*sqrt(-1120*I*z0 - 1680*I*z1 + 2800*I*z3 - 340*c0**2*z0**2 - 280*c0**2*z0*z1 + 40*c0**2*z0*z2 + 920*c0**2*z0*z3 + 420*c0**2*z1**2 + 120*c0**2*z1*z2 - 680*c0**2*z1*z3 + 12*c0**2*z2**2 - 184*c0**2*z2*z3 - 28*c0**2*z3**2 - 180*c0*c1*z0**2 - 240*c0*c1*z0*z1 + 600*c0*c1*z0*z3 + 180*c0*c1*z1**2 + 144*c0*c1*z1*z2 - 264*c0*c1*z1*z3 + 36*c0*c1*z2**2 - 216*c0*c1*z2*z3 - 60*c0*c1*z3**2 + 20*c0*c3*z0**2 + 72*c0*c3*z0*z1 + 88*c0*c3*z0*z2 - 200*c0*c3*z0*z3 - 12*c0*c3*z1**2 - 48*c0*c3*z1*z2 - 60*c0*c3*z2**2 + 80*c0*c3*z2*z3 + 60*c0*c3*z3**2 - 45*c1**2*z0**2 - 90*c1**2*z0*z1 - 18*c1**2*z0*z2 + 198*c1**2*z0*z3 + 27*c1**2*z1**2 + 54*c1**2*z1*z2 - 18*c1**2*z1*z3 + 27*c1**2*z2**2 - 90*c1**2*z2*z3 - 45*c1**2*z3**2 + 18*c1*c3*z0**2 + 84*c1*c3*z0*z1 + 120*c1*c3*z0*z2 - 240*c1*c3*z0*z3 + 18*c1*c3*z1**2 - 120*c1*c3*z1*z3 - 90*c1*c3*z2**2 + 60*c1*c3*z2*z3 + 150*c1*c3*z3**2 - 5*c3**2*z0**2 - 50*c3**2*z0*z1 - 230*c3**2*z0*z2 + 290*c3**2*z0*z3 - 45*c3**2*z1**2 - 270*c3**2*z1*z2 + 410*c3**2*z1*z3 + 75*c3**2*z2**2 + 350*c3**2*z2*z3 - 525*c3**2*z3**2))/(6*(2*z0 + 3*z1 - 5*z3))
        return solution.real

