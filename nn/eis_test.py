import numpy as np
import matplotlib.pyplot as plt
import os

from network import Network
from fc_layer import FCLayer
from activation_layer import ActivationLayer
from funcs import *

plt.style.use("C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle")
data_folder = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-04-19'


def extract_titration(d_folder):
    concs = []
    data = []
    
    for folder in os.listdir(d_folder):
        _, conc = folder.split('_')
        
        conc = float(conc)
        
        if conc != 0:
            for file in os.listdir(os.path.join(d_folder, folder)):
                if file.endswith('s.txt'):
                    file = os.path.join(d_folder, folder, file)
                    f, re, im = np.loadtxt(file, skiprows=1, unpack=True)
                    
                    # Sometimes 1Hz re is < 0 due to noise, throw out that data file
                    if not any([i<0 for i in re]):
                        concs.append(np.log10(conc)/10)
                        data.append([np.hstack((np.log10(re)/10, np.log10(-im)/10))])
    
                
    data  = np.array(data)
    concs = np.array(concs)
    
    return data, concs


def extract_invivo(d_folder):
    data  = []
    
    time_file = os.path.join(d_folder, '0000_time_list.txt')
    times = np.loadtxt(time_file)
    
    
    for file in os.listdir(d_folder):
        if file.endswith('s.txt'):
            file = os.path.join(d_folder, file)
            f, re, im = np.loadtxt(file, skiprows=1, unpack=True)
            
            data.append([np.hstack((np.log10(re)/10, np.log10(-im)/10))])
            
    
    return data, times



def train_network(d_folder, save_dir = None):
    # Generates a network and trains on titration data in d_folder
    
    data, concs = extract_titration(d_folder)

    net = Network()
    
    net.add(FCLayer(22*2, 100))
    net.add(ActivationLayer(tanh, tanh_prime))
    net.add(FCLayer(100, 80))
    net.add(ActivationLayer(tanh, tanh_prime))
    net.add(FCLayer(80, 50))
    net.add(ActivationLayer(tanh, tanh_prime))
    net.add(FCLayer(50, 1))
    net.add(ActivationLayer(tanh, tanh_prime))
    
    net.use(mse, mse_prime)
    net.fit(data, concs, epochs=10000, learning_rate=0.001, print_every=100)
    
    # Test network
    out = net.predict(data)
    vals = np.array([o[0][0] for o in out])
    plot_concs = concs
    
    # Plot test results
    fig, ax = plt.subplots(figsize=(5,5), dpi=100)
    ax.plot(10**(10*vals), label='Predicted')
    ax.plot(10**(10*plot_concs), label='Actual')
    ax.set_yscale('log')
    ax.set_ylabel('Concentration')
    ax.set_xlabel('File')
    ax.legend()
    
    # Plot error vs epoch number
    fig, ax = plt.subplots(figsize=(5,5), dpi=100)
    ax.plot(net.errs)
    ax.set_ylabel('log error')
    ax.set_xlabel('log epoch')
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    if save_dir:
        net.save(save_dir)
        print(f'Saved to {save_dir}')
    
    return net



def invivo(data_dir, net_dir = None, net = None):
    
    data, times = extract_invivo(data_dir)
    
    if net:
        net = net
    
    elif net_dir:
        net = Network(load_directory = net_dir)
        
    else:
        return
    
    concs = net.predict(data)
    concs = np.array([o[0][0] for o in concs])
    concs = 10**(10*concs)
    
    return concs, times
    

net = train_network(data_folder, save_dir=r'C:\Users\BRoehrich\Desktop\git\echem\nn\layers\phe_2')

concs, times = invivo(r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2022-05-11/rat 1 -330mV', 
                      net = net)

from scipy import signal

t = times[~np.isnan(concs)]
c = concs[~np.isnan(concs)]

f = signal.savgol_filter(c, 15, 1)

fig, ax = plt.subplots(figsize=(5,5), dpi=100)
ax.plot(t, c)
ax.plot(t, f)
