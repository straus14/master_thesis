__author__ = 'straus14'
import numpy as np
def db2pow(x):
    x_db = np.array(x)
    return (10**(x_db/10.0))

def dbm2pow(x):
    x_dbm = np.array(x)
    return ((10**(x_dbm/10.0))*(10**-3))




