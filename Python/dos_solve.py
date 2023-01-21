import numpy as np
from spetial_function import J_0,J_1,H1_1,H1_0
from make_param_vecs import make_param_vecs
from dos_quad import dos_quad
def Ldiag(t,tau,k, n, x, dx, ddx, f):
    si = np.where(t==tau)
    tsi = t[si]
    dxt = dx(tsi).T
    rotdxt = ddx(tsi).T@np.array(((0,-1),(1,0)))
    nrmdxt2 = np.sum(dxt**2,axis=1)
    v = 1/(2*np.pi)*np.sum(dxt*rotdxt,axis = 1)/nrmdxt2
    return v,si[0]
def L(t,tau,k, n, x, dx, ddx, f):
    xtau_t = x(tau)-x(t)
    xtau_t = xtau_t[:,:,0].T
    rotdxtau = dx(tau)[:,:,0].T@np.array(((0,-1),(1,0)))
    nrmxt_tau = np.sqrt(np.sum(xtau_t**2,axis = 1))
    v = 1j*k/2*np.sum(rotdxtau*xtau_t,axis=1)*H1_1(k*nrmxt_tau)/nrmxt_tau
    diagv,si = Ldiag(t,tau,k, n, x, dx, ddx, f)
    v[si] = diagv
    return v.reshape(len(v),1)
def L1(t,tau,k, n, x, dx, ddx, f):
    xt_tau = x(t)-x(tau)
    xt_tau = xt_tau[:,:,0].T
    rotdxtau = dx(tau)[:,:,0].T@np.array(((0,-1),(1,0)))
    nrmxt_tau = np.sqrt(np.sum(xt_tau**2,axis = 1))
    v = k/(2*np.pi)*np.sum(rotdxtau*xt_tau,axis = 1)*J_1(k*nrmxt_tau)/nrmxt_tau
    diagv,si = Ldiag(t,tau,k, n, x, dx, ddx, f)
    v[si] = np.zeros_like(diagv)
    return v.reshape(len(v),1)
def L2(t,tau,k, n, x, dx, ddx, f):
    v = L(t,tau,k, n, x, dx, ddx, f)-L1(t,tau,k, n, x, dx, ddx, f)*np.log(4*np.sin((t-tau)/2)**2)
    diagv,si = Ldiag(t,tau,k, n, x, dx, ddx, f)
    v[si,0] = diagv
    return v.reshape(len(v),1)
def Mdiag(t,tau,k, n, x, dx, ddx, f):
    si = np.where(t==tau)
    tsi = t[si]
    dxt = dx(tsi).T
    nrmdxt = np.sqrt(np.sum(dxt**2,axis=1))
    C = 0.57721
    v = (1j/2-C/np.pi-1/np.pi*np.log(k/2*nrmdxt))*nrmdxt
    return v,si[0]
def M(t,tau,k, n, x, dx, ddx, f):
    xt_tau = x(t)-x(tau)
    xt_tau = xt_tau[:,:,0].T
    nrmxt_tau = np.sqrt(np.sum(xt_tau**2,axis=1))
    dxtau = dx(tau)[:,:,0].T
    nrmdxtau = np.sqrt(np.sum(dxtau**2,axis = 1))
    v = 1j/2*H1_0(k*nrmxt_tau)*nrmdxtau
    return v.reshape(len(v),1)
def M1(t,tau,k, n, x, dx, ddx, f):
    xt_tau = x(t)-x(tau)
    xt_tau = xt_tau[:,:,0].T
    nrmxt_tau = np.sqrt(np.sum(xt_tau**2,axis=1))
    dxtau = dx(tau)[:,:,0].T
    nrmdxtau = np.sqrt(np.sum(dxtau**2,axis=1))
    v = -1/(2*np.pi)*J_0(k*nrmxt_tau)*nrmdxtau
    return v.reshape(len(v),1)
def M2(t,tau,k, n, x, dx, ddx, f):
    v = M(t,tau,k, n, x, dx, ddx, f)-M1(t,tau,k, n, x, dx, ddx, f)*np.log(4*np.sin((t-tau)/2)**2)
    diagv,si = Mdiag(t,tau,k, n, x, dx, ddx, f)
    v[si,0] = diagv
    return v
def getK(t,k, n, x, dx, ddx, f,eta):
    tau = np.linspace(0,2*n-1,2*n)*np.pi/n
    t,tau,s = make_param_vecs(t,tau)
    K1 = L1(t,tau,k, n, x, dx, ddx, f)+1j*eta*M1(t,tau,k, n, x, dx, ddx, f)
    K2 = L2(t,tau,k, n, x, dx, ddx, f)+1j*eta*M2(t,tau,k, n, x, dx, ddx, f)
    K1 = K1.reshape(s)
    K2 = K2.reshape(s)
    return K1,K2
def getg(t,k, n, x, dx, ddx, f):
    g = np.zeros((len(t),1))+1j*np.zeros((len(t),1))
    for ii in range(len(t)):
        g[ii] = 2*f(t[ii])
    return g
def setup_sovler(k, n, x, dx, ddx, f,eta):
    t = np.linspace(0,2*n-1,2*n)*np.pi/n
    K1,K2 = getK(t,k, n, x, dx, ddx, f,eta)
    R = dos_quad(n,t)
    A = R*K1+np.pi/n*K2
    g = getg(t,k, n, x, dx, ddx, f)
    OP = np.eye(A.shape[0])-A
    return OP,t,g,R,A
def solve_for_phi(k, n, x, dx, ddx, f,eta):
    OP,t,g,R,A = setup_sovler(k, n, x, dx, ddx, f,eta)
    phi = np.linalg.solve(OP,g)
    return phi,OP,t,g,R,A