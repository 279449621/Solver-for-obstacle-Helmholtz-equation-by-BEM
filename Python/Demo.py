import numpy as np
from dos_solve import solve_for_phi
from find_u import find_u
from find_far_field import find_far_field
x = lambda t:np.array((-0.65+np.cos(t)+0.65*np.cos(2*t),1.5*np.sin(t)))
dx = lambda t:np.array((-np.sin(t)-1.3*np.sin(2*t),1.5*np.cos(t)))
ddx = lambda t:np.array((-np.cos(t)-2.6*np.cos(2*t),-1.5*np.sin(t)))
k = 1
inc_ang = 0
inc_dir = np.array((np.cos(np.pi*inc_ang/180),np.sin(np.pi*inc_ang/180)))
inc_dir.reshape((2,1))
f = lambda t:-np.exp(1j*k*x(t)@inc_dir)
n = 32
eta = k
#Solve Helemhotz equation for phi
phi,OP,t,g,R,A = solve_for_phi(k, n, x, dx, ddx, f,eta)

#Solve Helemhotz equation for u
xn = 200
yn = 200
xs = np.linspace(-8,8,xn)
ys = np.linspace(-8,8,yn)
ts = t
v,xv,yv = find_u(xs,ys,ts,k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A)

#Find u_inf
angles = np.linspace(0,2*np.pi,129)
angles = angles[0:-1]
uinf,nu,y,xhat,phi = find_far_field(k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A,angles)
print(uinf)
