import numpy as np

# Activation functions
def tanh(x):
    return np.tanh(x)


def tanh_prime(x):
    return 1-np.tanh(x)**2


def elu(z,alpha=-1):    
	return np.where(z>0, z, np.exp(z)-1)

def elu_prime(z,alpha=-1):
	return np.where(z>0, 1, np.exp(z))







# Loss function
def mse(y_true, y_pred):
    return np.mean(np.power(y_true-y_pred, 2))

# Loss function derivative
def mse_prime(y_true, y_pred):
    return 2*(y_pred-y_true)/y_true.size
