from layer import Layer


class ActivationLayer(Layer):
    
    def __init__(self, activation, activation_prime):
        self.activation       = activation
        self.activation_prime = activation_prime
        
        
    def forward_propagation(self, input_data):
        # return activated input f(X)
        self.input  = input_data
        self.output = self.activation(self.input)
        return self.output
    
    
    def backward_propagation(self, output_error, learning_rate):
        # Returns input_error = dE/dX for given output_error=dE/dY
        return self.activation_prime(self.input) * output_error