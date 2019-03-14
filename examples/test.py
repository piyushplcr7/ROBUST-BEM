import numpy as np
from numpy import linalg as LA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from math import log

class lineparam:
    def __init__(self,x1,y1,x2,y2):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

    def point(self,t):
        meanx = (self.x2+self.x1)/2.
        meany = (self.y2+self.y1)/2.
        diffx = (self.x2-self.x1)/2.
        diffy = (self.y2-self.y1)/2.
        return np.array([t*diffx+meanx,t*diffy+meany])

class circleparam:
    def __init__(self,x1,y1,R,phi_s,phi_e):
        self.x1 = x1
        self.y1 = y1
        self.r = R
        self.phi_s = phi_s
        self.phi_e = phi_e

    def point(self,t):
        meanphi = (self.phi_s+self.phi_e)/2.
        difffi = (self.phi_e-self.phi_s)/2.
        phi = t*difffi + meanphi
        r = self.r
        x = self.x1
        y = self.y1
        return np.array([x+r*np.cos(phi),y+r*np.sin(phi)])

class paraparam:
    def __init__(self,a):
        self.a = a

    def point(self,t):
        return np.array([t**2+self.a,t])

def func(curve1,s,curve2,t,eps):
    blowup = 1/((s-t)**2+eps**2)#1/(eps**2) - (1/eps**2-1) * (s-t)**2 * 0.25
    #blowup = 1/((abs(s-t)+eps)**2)
    distance2 = LA.norm(curve1.point(s)-curve2.point(t))**2
    return blowup * distance2
    #return distance2


eps = .01
#curve1 = lineparam(0,0,-eps,1)
#curve2 = lineparam(eps,0,2*eps,1)
#curve1 = circleparam(0,0,1,-np.pi/2,np.pi/2)
#curve2 = circleparam(eps,0,1,-np.pi/2,np.pi/2)
curve1 = paraparam(0)
curve2 = paraparam(eps)
N = 100
s = np.linspace(-1,1,N)
t = np.linspace(-1,1,N)
x,y = np.meshgrid(s,t)
z = np.empty(x.shape)
for i in np.arange(len(s)):
    for j in np.arange(len(t)):
        z[i,j] = func(curve1,s[i],curve2,s[j],eps)

print ("minimum: ", np.amin(z))
fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel('s')
ax.set_ylabel('t')

plt.figure()
plt.plot(s,z[49,:])
plt.show()
