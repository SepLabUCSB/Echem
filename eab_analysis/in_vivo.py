import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from EIS_Fit import EIS_fit

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# data_dir = r'C:/Users/BRoehrich/Desktop/2022-04-13/2022-04-13/6 hr with dosing -330mV 1154am'
data_dir = r'C:/Users/BRoehrich/Desktop/2022-04-13/2022-04-13/first 200 fit testing'
# data_dir = r'C:\Users\BRoehrich\Desktop\2022-04-13\2022-04-13\PBS 0mM target'
# data_dir = r'C:\Users\BRoehrich\Desktop\2022-04-13\2022-04-13\30 min equilibration -350mV 1120am'


target = 'Phenylalanine'

starting_bounds = {
    'R1': [1e-1, 1e3],
    'R2': [1000, 20000],
    'Q1': [1e-10, 1e-5],
    'n1': [1,1],
    'Q2': [1e-15, 1e-5],
    'n2': [0.9,1.1]
    }

absolute_bounds = {
    'R1': [1e-2, 1e4],
    'R2': [1000, 75000],
    'Q1': [1e-10, 1e-4],
    'n1': [1,1],
    'Q2': [1e-15, 1e-4],
    'n2': [0.9,1.1]
    }

# Vancomycin
# starting_guess = {
#     'R1': 284, 
#     'R2': 40000, 
#     'Q1': 1.8e-07, 
#     'n1': 1, 
#     'Q2': 3.2e-07, 
#     'n2': 0.9
#     }


# Tobramycin
# starting_guess = {
#     'R1': 100, 
#     'R2': 10000, 
#     'Q1': 5e-07, 
#     'n1': 1, 
#     'Q2': 5e-07, 
#     'n2': 1
#     }


# Phenylalanine
starting_guess = {
    'R1': 130, 
    'R2': 7000, 
    'Q1': 6e-07, 
    'n1': 1, 
    'Q2': 1e-06, 
    'n2': 0.95
    }


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
        
        self.params = params
        self.circuit = 'Randles_adsorption'
        
    
    def fit(self, starting_guess = None,
            plot = False, bounds=starting_bounds,
            ga = True, LEVM = True, **kwargs):
        
        DataFile = EIS_fit.DataFile(self.file, circuit=self.circuit, 
                                    Z=self.Z, freqs = self.freqs,
                                    bounds=bounds, params = self.params)
        
        if ga:
            DataFile.ga_fit(starting_guess = starting_guess, **kwargs)
        
        else:
            DataFile.params = starting_guess
           
        
        if LEVM:
            DataFile.LEVM_fit()
        
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
        
        

def compress_data(output_folder, freqs, re, im, num):   
    d = pd.DataFrame(
            {'f': freqs,
            're': re,
            'im': im}
            )
    
    fname = output_folder + f'\\{num}s.txt'
    
    d.to_csv(fname, columns = ['f', 're', 'im'],
                 header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], 
                 sep = '\t', index = False, encoding='ascii')
        
 
    
def extract_data(folder, average = 1, freq_cutoff = 100000):
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
            df = df[df['f'] <= freq_cutoff]
            
            res.append(df['re'].to_numpy())
            ims.append(df['im'].to_numpy())
            fs.append(df['f'].to_numpy())
            
            
            
            
            
    output_folder = 'C:/Users/BRoehrich/Desktop/compressed'      
    
    specs = []
    # Average over n spectra
    for i in range(int(len(files)/average)):
        x = i*average
        y = (i+1)*average
        spec = Spectrum(ts[x], fs[x], 
                        np.mean(res[x:y], axis=0), 
                        np.mean(ims[x:y], axis=0),
                        files[x])
        specs.append(spec)
        
        # compress_data(output_folder, fs[x], 
        #                 np.mean(res[x:y], axis=0), 
        #                 np.mean(ims[x:y], axis=0),
        #                 files[x][:4])
        
    
    # Running average over previous n spectra
    # for i in range(average, len(files)):
    #     x = i - average
        
    #     spec = Spectrum(ts[i], fs[i],
    #                     np.mean(res[x:i], axis=0),
    #                     np.mean(ims[x:i], axis=0),
    #                     files[i])
    #     specs.append(spec)
        
    
    print('Loaded all files.')
    
    return specs



def fit_all_data(specs):
    # specs: list of Spectrum objects
    
    specs.sort(key=lambda s:s.time)
    
    for i in range(len(specs)):
        
        if i == 0:
            specs[i].fit(ga=False, LEVM=True, starting_guess = starting_guess, n_iter = 50)
            continue
            
        if i < 100:
            specs[i].fit(ga=False, LEVM=True, starting_guess = specs[i-1].params)
            continue
        
        # plot = False
        
        # if i == 0:
        #     n_iter = 200
        #     r_mut = 0.4
        #     bounds = starting_bounds
        #     guess = starting_guess
        #     plot = True
        #     ga = True
        #     LEVM = True
        
        # else:
        #     n_iter = 50
        #     r_mut = 0.25
        #     ga = False
        #     LEVM = True
        #     guess = specs[i-1].params
        #     # if i %100 == 0:
        #     #     plot = True
            
        #     bounds = {param:[0.8*val, 1.2*val] 
        #                  for param, val in specs[i-1].params.items()}
            
        #     for param, val in bounds.items():
        #         if param.startswith('n'):
        #             bounds[param] = absolute_bounds[param]
        #         if val[0] < absolute_bounds[param][0]:
        #             bounds[param][0] = absolute_bounds[param][0]
        #         if val[1] > absolute_bounds[param][1]:
        #             bounds[param][1] = absolute_bounds[param][1]
                
        
        
        # specs[i].fit(n_iter = n_iter, starting_guess = guess, 
        #              plot = plot, bounds = bounds, 
        #              ga = ga, LEVM = LEVM, r_mut = r_mut)
            

            
        # if specs[i].chi_squared > 90:     
        #     print('Poor fit, restarting')
        #     specs[i].fit(n_iter = 50, starting_guess = specs[i-2].params, plot=True,
        #               bounds = absolute_bounds, ga=True, LEVM=True)
        #     specs[i].params = specs[i-1].params
        #     specs[i].chi_squared = specs[i-1].chi_squared
        # if not any([absolute_bounds[param][0] < 
        #         specs[i].params[param] < 
        #         absolute_bounds[param][1]
        #         for param in specs[i].params]):
        #     specs[i].params = specs[i-1].params
        #     specs[i].chi_squared = specs[i-1].chi_squared
            
        #     print([param for param in specs[i].params if absolute_bounds[param][0] < 
        #         specs[i].params[param] < 
        #         absolute_bounds[param][1]])
            
        #     print('Bad parameter value, restarting')
        #     specs[i].fit(n_iter = 50, starting_guess = specs[i-2].params, plot=True,
        #               bounds = absolute_bounds, ga=True, LEVM=True)

















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
    


def thresholding_algo(y, lag, threshold, influence):   
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
    
    
    return signals, np.array(avgFilter), np.array(stdFilter)



def plot_filtered_ket(specs):
    ts = []
    kets = []

    for spec in specs:
        ts.append(spec.time)
        kets.append(spec.ket)
    
    # y = signal.savgol_filter(kets, 15, 1)
    signals, y, std = thresholding_algo(kets, 5, 1, 0.3)
    
    fig, ax = plt.subplots()
    ax.plot(ts, kets)
    ax.plot(ts, y)
    # # Vertical lines to show injection points
    # ax.axvline(1640, 0.1, 0.8, color='k', linestyle='--')
    # ax.axvline(8200, 0.1, 0.8, color='k', linestyle='--')
    ax.set_title('ket')
    ax.set_xlabel('Time/ s')
    ax.set_ylabel('$k_{et}/  s^{-1}$')



def plot_freq_vs_t(specs, freq):
    
    i = np.where(specs[0].freqs.astype(int) == int(freq))[0][0]
    
    ts = []
    vals = []
    
    for spec in specs:
        ts.append(spec.time)
        vals.append(spec.phase[i])
    
    fig, ax = plt.subplots()
    ax.plot(ts, vals)
    ax.set_title(f'Phase @ {int(freq)} Hz')
    ax.set_xlabel('Time/ s')
    ax.set_ylabel('Phase/ $\degree$')
    



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
    
specs = extract_data(data_dir, average=1, freq_cutoff=30000)
# plot_phase_maps(specs)
fit_all_data(specs)
# plot_params(specs)
# plot_filtered_ket(specs)




    
    
    
    
    
    
    
    
    
    
    