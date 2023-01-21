import numpy as np
def make_param_vecs(t,tau):
    t = t.reshape(len(t),1)
    tau = tau.reshape(1,len(tau))
    s = t.shape
    t = np.tile(t,tau.shape)
    tau = np.tile(tau,s)
    matsize = t.shape
    vecsize = int(matsize[0]*matsize[1])
    tau = tau.reshape(vecsize,1)
    t = t.reshape(vecsize,1)
    return t,tau,matsize