import numpy as np

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def polar_to_cartesian(z,theta):
    return z*np.exp(1j*np.deg2rad(theta))

def cartesian_to_polar(Z):
    x = Z.real
    y = Z.imag
    z = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return z,theta

def MAPE(predictions, targets): 
    targets, predictions = np.array(targets), np.array(predictions)
    return np.mean(np.abs((targets - predictions) / targets)) * 100