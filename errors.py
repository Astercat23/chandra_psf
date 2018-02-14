import numpy as np

def prop_add(e1, e2):
    return np.sqrt(e1**2 + e2**2)

def prop_mult(x, y, x_err, y_err):
    # Error propogation for F = x * y
    xy = x * y
    return xy * np.sqrt((x_err/x)**2 + (y_err/y)**2)

def prop_div(x, y, x_err, y_err):
    # Error propogation for F = x / y
    xy_div = x / y
    return xy_div * np.sqrt((x_err/x)**2 + (y_err/y)**2)
