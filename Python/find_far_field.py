import numpy as np
def unit_normal(k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A,angles):
    v = dx(t)[:,:,0].T@np.array(((0,-1),(1,0)))
    r = np.sqrt(np.sum(v**2,axis=1))
    r = r.reshape(len(r),1)
    v = v/np.tile(r,(1,2))
    return v
def find_far_field(k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A,angles):
    angles = angles.reshape(1,len(angles))
    t = t.reshape(len(t),1)
    y = x(t)[:,:,0].T
    dy = dx(t)[:,:,0].T
    ds = np.sqrt(np.sum(dy**2,axis=1))
    phi = phi.reshape(phi.shape[0])
    xhat = np.concatenate((np.cos(angles),np.sin(angles)),axis=0)
    nu = unit_normal(k, n, x, dx, ddx, f,eta,phi,OP,t,g,R,A,angles)
    uinf = np.zeros((angles.shape[1],1))+1j*np.zeros((angles.shape[1],1))
    a = np.exp(-1j*np.pi/4)/np.sqrt(8*np.pi*k)
    for j in range(xhat.shape[1]):
        xh = xhat[:,j]
        K = (k*nu@xh+eta)*np.exp(-1j*k*y@xh)*phi*ds
        uinf[j] = a*np.sum(K*np.pi/n,axis=0)
    return uinf,nu,y,xhat,phi