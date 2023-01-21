import numpy as np
def dos_quad(n,tn):
    tn = tn.reshape(len(tn),1)
    tj = np.linspace(0,2*n-1,2*n)*np.pi/n
    tj = tj.reshape(len(tj),1)
    m = np.linspace(1,n-1,n-1)
    m = m.reshape(1,len(m))
    R = np.zeros((tj.size,tj.size))
    for ii in range(len(tn)):
        ti = tn[ii]
        M = np.cos(m*(ti-tj))/np.tile(m,np.shape(tj))
        M = M.T
        temp = np.pi/n**2*np.cos(n*(ti-tj))
        temp = temp.reshape(len(temp))
        R[ii,:] = -2*np.pi/n*np.sum(M,axis = 0)-temp
    return R