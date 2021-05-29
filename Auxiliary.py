import numpy as np

def u_m(m,full_supply,eta_max_m,eta_0):
    if(full_supply):
        fac = np.ones_like(m);
    else:
        fac = np.tanh(m/eta_max_m*eta_0);

    return fac*eta_max_m;

def u_p(p,full_supply,eta_max_p,eta_0):

    if(full_supply):
        fac = np.ones_like(p);
    else:
        fac = np.tanh(p/eta_max_p*eta_0);

    return fac*eta_max_p;
