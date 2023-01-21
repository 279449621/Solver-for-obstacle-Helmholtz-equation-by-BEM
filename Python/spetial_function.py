from scipy.special import jv,hankel1
def J_0(x):
    return jv(0,x)
def J_1(x):
    return jv(1,x)
def H1_0(x):
    return hankel1(0,x)
def H1_1(x):
    return hankel1(1,x)