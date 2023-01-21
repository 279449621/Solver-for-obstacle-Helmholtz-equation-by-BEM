import numpy as np
from spetial_function import H1_0,H1_1,J_0,J_1
from make_param_vecs import make_param_vecs
def Loffbdy(xv,tau,k, n, x, dx, ddx, f,eta):
    xt_tau = x(tau)[:,:,0].T-xv
    rotdxtau = dx(tau)[:,:,0].T@np.array(((0,-1),(1,0)))
    nrmxt_tau = np.sqrt(np.sum(xt_tau**2,axis = 1))
    v = -1j*k/4*np.sum(rotdxtau*xt_tau,axis = 1)*H1_1(k*nrmxt_tau)/nrmxt_tau
    return v
def Moffby(xv,tau,k, n, x, dx, ddx, f,eta):
    xt_tau = xv-x(tau)[:,:,0].T
    nrmxt_tau = np.sqrt(np.sum(xt_tau**2,axis=1))
    dxtau = dx(tau)[:,:,0].T
    nrmdxtau = np.sqrt(np.sum(dxtau**2,axis=1))
    v = -1j/4*H1_0(k*nrmxt_tau)*nrmdxtau
    return v
def find_u(xs,ys,ts,k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A):
    xv,yv,temp = make_param_vecs(xs,ys)
    xvv,temp,s1 = make_param_vecs(xv,ts)
    yvv,tauv,temp = make_param_vecs(yv,ts)
    vn = np.concatenate((xvv,yvv),axis=1)
    Koffbdy = (Loffbdy(vn,tauv,k, n, x, dx, ddx, f,eta)+1j*eta*Moffby(vn,tauv,k, n, x, dx, ddx, f,eta))*np.pi/n
    Koffbdy = Koffbdy.reshape(s1)
    v = Koffbdy@phi
    return v,xv,yv
