import numpy as np
import os

from fc_layer import FCLayer
from activation_layer import ActivationLayer
from funcs import tanh, tanh_prime, mse, mse_prime, elu, elu_prime


func_dict = {
    'tanh': tanh,
    'tanh_prime': tanh_prime,
    'elu': elu,
    'elu_prime': elu_prime,
    'mse': mse,
    'mse_prime': mse_prime 
    }



class Network:
    
    def __init__(self, load_directory = None):
        self.layers = []
        self.loss = None
        self.loss_prime = None
        self.errs = []
        
        if load_directory:
            self.load_from_dir(load_directory)
    

        
    # Add layer to network
    def add(self, layer):
        self.layers.append(layer)
        
        
    # set loss to use
    def use(self, loss, loss_prime):
        self.loss = loss
        self.loss_prime = loss_prime
       
        
    # predict output for given input
    def predict(self, input_data):
        result = []
        
        for i in range(len(input_data)):
            output = input_data[i]
            for layer in self.layers:
                output = layer.forward_propagation(output)
            result.append(output)
            
        return result
    
    
    # Train network
    def fit(self, x_train, y_train, epochs, learning_rate, print_every=1):
        samples = len(x_train)
        
        
        for i in range(epochs):
            err = 0
            for j in range(samples):
                # Forward propagate
                output = x_train[j]
                for layer in self.layers:
                    output = layer.forward_propagation(output)
                
                
                # Calculate loss using knowns (y_train)
                err += self.loss(y_train[j], output)
                
                # Backwards propagation
                error = self.loss_prime(y_train[j], output)
                
                for layer in reversed(self.layers):
                    error = layer.backward_propagation(error, learning_rate)
                
            # Get average error for this epoch
            err = err/samples
            self.errs.append(err)
            if i%print_every == 0:
                print(f'Epoch {i}/{epochs}: error={err}')
                
    
    def save(self, directory):
        
        txtfile = os.path.join(directory, 'params.txt')
        with open(txtfile, 'w') as f:
            for layer in self.layers:
                if hasattr(layer, 'activation'):
                    f.write(str(layer.activation).split(' ')[1])
                    f.write(',')
                    f.write(str(layer.activation_prime).split(' ')[1])
                    f.write('\n')
            f.write(str(self.loss).split(' ')[1])
            f.write(',')
            f.write(str(self.loss_prime).split(' ')[1])
        
        
        i = 0
        for layer in self.layers:
            if hasattr(layer, 'weights'):
                fname = os.path.join(directory, f'{i}.csv')
                arr = np.vstack((layer.weights, layer.bias))
                np.savetxt(fname, arr, delimiter=',')
                
                i += 1
   
                
    def load_from_dir(self, directory):
        # Initialize pretrained network
        
        # Get activation and loss functions
        txtfile = os.path.join(directory, 'params.txt')
        l = []
        with open(txtfile, 'r') as f:
            for line in f:
                l.append(tuple(line.strip('\n').split(',')))
                
        # Get layer weight and bias matrices
        layerfiles = [os.path.join(directory, file)
                      for file in os.listdir(directory)
                      if file.endswith('.csv')]
        
        i = 0
        for file in layerfiles:
            arr = np.loadtxt(file, delimiter=',')
            weights = arr[:-1]
            bias = arr[-1:]
            input_size = weights.shape[0]
            try:
                output_size = weights.shape[1]
            except:
                output_size = 1
            # Add FCLayer with this matrix
            self.add(FCLayer(input_size, output_size, weights, bias))
            
            # Corresponding activation functions
            activations = (
                            func_dict[l[i][0]],
                            func_dict[l[i][1]]
                           )
            
            
            self.add(ActivationLayer(*activations))
            
            i += 1
        
        loss_funcs = (
                     func_dict[l[i][0]],
                     func_dict[l[i][1]]
                    )
        
        self.use(*loss_funcs)
                
            
            
            
            
            
            
            
            