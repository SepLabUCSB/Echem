import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from EIS_Fit import EIS_fit

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


data_dir = r'C:/Users/BRoehrich/Desktop/2022-04-05/2022-04-05/long time dosing'

target = 'Vancomycin'

starting_bounds = {
    'R1': [1e-1, 1e9],
    'R2': [10000, 150000],
    'Q1': [1e-15, 1],
    'n1': [0.8,1.1],
    'Q2': [1e-15, 1],
    'n2': [0.5,1.2]
    }

# Restrictive Rct bounds
# bounds = {
#     'R1': [1e-1, 1e9],
#     'R2': [1, 50000],
#     'Q1': [1e-15, 1],
#     'n1': [0.8,1.1],
#     'Q2': [1e-15, 1],
#     'n2': [0.8,1.1]
#     }

# Vancomycin
starting_guess = {
    'R1': 284, 
    'R2': 40000, 
    'Q1': 1.8e-07, 
    'n1': 1, 
    'Q2': 3.2e-07, 
    'n2': 0.9
    }


# Tobramycin
# starting_guess = {
#     'R1': 100, 
#     'R2': 10000, 
#     'Q1': 5e-07, 
#     'n1': 1, 
#     'Q2': 5e-07, 
#     'n2': 1
#     }

units = {
    'R1': '$\Omega$',
    'R2': '$\Omega$',
    'Q1': 'F',
    'Q2': 'F',
    'n1': '',
    'n2': ''    
    }


class Spectrum:
        
    def __init__(self, time, freqs, re, im, file, AC = None,
                 DC = None, elec_size = None, params = None):
        
        self.file   = file
        self.time   = float(time)
        self.freqs  = freqs.astype(float)
        self.re     = re
        self.im     = im
        self.AC     = AC
        self.DC     = DC
        self.elec_size = elec_size
        
        self.Z = re + 1j*im
        self.phase = np.angle(self.Z, deg=True)
        
        # phase_peak_i = np.where(self.phase == max(self.phase[:-4]))[0][0]
        # phase_min_i = np.where(self.phase == min(self.phase[6:]))[0][0]
        # self.phase_peak = self.freqs[phase_peak_i]
        # self.phase_min = self.freqs[phase_min_i]
        
        self.params = params
        
    
    def fit(self, n_iter=None, starting_guess = None,
            plot = False, bounds=starting_bounds, **kwargs):
        self.circuit = 'Randles_adsorption'
        DataFile = EIS_fit.DataFile(self.file, circuit=self.circuit, 
                                    Z=self.Z, freqs = self.freqs,
                                    bounds=bounds)
        
        
        DataFile.ga_fit(n_iter = n_iter, starting_guess = starting_guess, **kwargs)
        t = DataFile.params
        DataFile.LEVM_fit()
        if DataFile.params == t:
            print('LEVM fit returns same params')
        
        if plot:
            fig, ax = plt.subplots()
            DataFile.plot_fit(ax=ax, Bode=True)
            plt.title(self.file)
            plt.show()
        
        # Get fit parameters
        self.params = DataFile.params
        self.ket = 1/(2*self.params['R2']*self.params['Q2'])
        self.chi_squared = DataFile.chi_squared
        print(self.params, self.chi_squared)
        
        
        
def extract_data(folder, average = 1):
    # Extract time series of impedance spectra
    
    res = []
    ims = []
    fs  = []
    files = []
    
    
    for file in os.listdir(folder):
        if file.endswith('time_list.txt'):
            file = os.path.join(folder, file)
            with open(file, 'r') as f:
                ts = f.read().split('\n')
        if file.endswith('fits.txt'): 
            continue
        if file.endswith('s.txt'):
            files.append(file)
            file = os.path.join(folder, file)
            df = pd.read_csv(file, skiprows=1, names=('f', 're', 'im'), sep='\t')
            res.append(df['re'].to_numpy())
            ims.append(df['im'].to_numpy())
            fs.append(df['f'].to_numpy())
            
    
    specs = []
    for i in range(int(len(files)/average)):
        x = i*average
        y = (i+1)*average
        try:
            spec = Spectrum(ts[x], fs[x], 
                            np.mean(res[x:y], axis=0), 
                            np.mean(ims[x:y], axis=0),
                            files[x])
            specs.append(spec)
        except:
            print(i)
            print(res[i])
            return
    
    return specs



def fit_all_data(specs):
    # specs: list of Spectrum objects
    
    specs.sort(key=lambda s:s.time)
    
    for i in range(len(specs)):
        if i == 0:
            specs[i].fit(n_iter = 200, starting_guess = starting_guess, plot=True,
                     bounds = starting_bounds)
            
        elif i%100 == 0:
            
            bounds = {param:[val/2, 2*val] 
                         for param, val in specs[i-1].params.items()}
            bounds['n1'] = [0.8,1.1]
            bounds['n2'] = [0.5,1.2]
            
            specs[i].fit(n_iter = 20, starting_guess = specs[i-1].params, plot=True,
                     bounds = bounds)
        
        else:
            bounds = {param:[val/2, 2*val] 
                         for param, val in specs[i-1].params.items()}
            bounds['n1'] = [0.8,1.1]
            bounds['n2'] = [0.5,1.2]
            
            
            specs[i].fit(n_iter = 20, starting_guess = specs[i-1].params, plot=False,
                     bounds = bounds)
            
        if specs[i].chi_squared > 25:
            bounds = {param:[val/2, 2*val] 
                         for param, val in specs[i-1].params.items()}
            bounds['n1'] = [0.8,1.1]
            bounds['n2'] = [0.5,1.2]
            
            specs[i].fit(n_iter = 200, starting_guess = starting_guess, plot=True,
                     bounds = bounds)



def plot_phase_maps(specs):
    colors = plt.cm.viridis(np.linspace(0.2,0.8,len(specs)))
    fig, ax = plt.subplots()
    for i in range(len(specs)):
        ax.plot(specs[i].freqs, specs[i].phase, color=colors[i])
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')
    ax.set_xscale('log')
    ax.set_xticks([1,10,1000])
            
    


def plot_params(specs):
    ts = []
    kets = []

    for spec in specs:
        ts.append(spec.time)
        kets.append(spec.ket)
    
    for param, _ in specs[0].params.items():
        l = []
        for spec in specs:
            l.append(spec.params[param])
        
        fig, ax = plt.subplots()
        ax.plot(ts, l)
        ax.set_title(param)
        ax.set_xlabel('Time/ s')
        ax.set_ylabel(f'{param}/ {units[param]}')
    
    fig, ax = plt.subplots()
    ax.plot(ts, kets)
    ax.set_title('ket')
    ax.set_xlabel('Time/ s')
    ax.set_ylabel('$k_{et}/  s^{-1}$')
    


def plot_filtered_ket(specs):
    ts = []
    kets = []

    for spec in specs:
        ts.append(spec.time)
        kets.append(spec.ket)
    
    y = signal.savgol_filter(kets, 15, 3)
    
    fig, ax = plt.subplots()
    ax.plot(ts, kets)
    ax.plot(ts, y)
    ax.set_title('ket')
    ax.set_xlabel('Time/ s')
    ax.set_ylabel('$k_{et}/  s^{-1}$')



def save_fits(specs):
    '''
    specs = list of Spectrum objects
    '''
        
    rows_list = []
    
    for spectrum in specs:
        dict1 = {
               'time': spectrum.time,
               'R1': spectrum.params['R1'],
               'R2': spectrum.params['R2'],
               'Q1': spectrum.params['Q1'],
               'n1': spectrum.params['n1'],
               'Q2': spectrum.params['Q2'],
               'n2': spectrum.params['n2'],
               'ket': spectrum.ket,
               'chi_squared': spectrum.chi_squared
               }
        
        rows_list.append(dict1)
        rows_list.sort(key= lambda d:d['time'])
        
        # out_df = out_df.append(data_list, ignore_index=True)
        df = pd.DataFrame(rows_list)
    
    df.to_csv(r'C:\Users\BRoehrich\Desktop\saved_fits.csv')
    
    return df    
    
# specs = extract_data(data_dir, average=5)
# plot_phase_maps(specs)
# fit_all_data(specs)
# plot_params(specs)




    
    
    
    
    
    
    
    
    
    
    