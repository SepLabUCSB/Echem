import numpy as np
from EIS_Fit import ga, circuits, LEVM
# import EIS_fit.circuits as circuits
# import EIS_fit.LEVM as LEVM
import matplotlib.pyplot as plt
import pandas as pd
import os
import time
# plt.style.use('Z:/Projects/Brian/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

'''
Based on ChemElectroChem 2021, 8, 2956–2967

doi.org/10.1002/celc.202100778
'''


path = r'C:/Users/BRoehrich/Desktop/EIS fit folder/2021-06-01_FcMeOH_50 mV_beads_009 (2)'

# circuit = 'Randles_uelec'
# bounds = {
#     'R1': [1e-1, 1e9],
#     'R2': [1e-1, 1e9],
#     'R3': [1e-1, 1e9],
#     'Q1': [1e-15, 1],
#     'n1': [1,1],
#     'Q2': [1e-15, 1],
#     'n2': [0,1]
#     }



'''
User inputs:
    circuit string
    bounds
    tolerence for fit

    folder containing all files to fit


Read data from file: freqs, re, im

Do genetic algorithm: fitted parameters

Do least squares fitting: refined parameters

Check if least squares fitting satisfies fit tolerence

No:
    Repeat genetic algorithm


Go to next file in batch, take refined parameters as initial values
for least squares fit


End: export csv of frequency, param 1, param 2, ...
            file 1
            file 2
            file 3
            ...

'''



class DataFile:
    
    def __init__(self, file, circuit, bounds, ax=None, Z=None,
                 freqs=None):
        
        if len(Z) == 0 and file != '':
            df = pd.read_csv(file, skiprows=1, names=
                           ('freqs', 're', 'im'), sep='\t')
            self.freqs = df['freqs'].to_numpy()
            self.re = df['re'].to_numpy()
            self.im = df['im'].to_numpy()
            self.Z = self.re + 1j*self.im
        
        else:
            self.Z = Z
            self.freqs = freqs
            self.re = np.real(Z)
            self.im = np.imag(Z)
        
        self.file = file
        self.circuit = circuit
        self.bounds = bounds
                        
        self.params = dict()
        self.score = 1e9
        
        self.ax = ax
            
    
        
    def ga_fit(self, ax=None, starting_guess=None, **kwargs):
        '''
        Perform genetic algorithm fit using ga.py
        '''
        if ax is not None:
            ax.plot(self.re/1e6, -self.im/1e6, 'o')
            
        self.params, self.score = ga.genetic_algorithm(
            self.freqs, self.Z, self.bounds, self.circuit, ax=ax,
            starting_guess = starting_guess, **kwargs)
        
        self.chi_squared = circuits.calc_chi_squared(self.freqs, self.Z, 
                                                     self.params, self.circuit)
        
        circuit_list, circuit_funcs = circuits.list_circuits()
        circuitfunc = circuit_funcs[circuit_list.index(self.circuit)]
        
        self.fits = circuitfunc(self.freqs, self.params)
        
        
    def LEVM_fit(self, timeout=2, **kwargs):
        '''
        Perform least-squares fit using LEVM.py
        '''
        try:
            self.params = LEVM.LEVM_fit(self.freqs, self.Z, 
                                        self.params, self.circuit,
                                        timeout = timeout)
        
        except:
            print('LEVM fit timed out, performing GA fit. File: ', 
                  self.file)
            self.ga_fit(starting_guess=self.params, n_iter = 50)
            pass
                
        self.score = circuits.leastsq_errorfunc(self.freqs, self.Z,
                                                self.params, self.circuit)
        
        self.chi_squared = circuits.calc_chi_squared(self.freqs, self.Z, 
                                                     self.params, self.circuit)
        
        circuit_list, circuit_funcs = circuits.list_circuits()
        circuitfunc = circuit_funcs[circuit_list.index(self.circuit)]
        
        self.fits = circuitfunc(self.freqs, self.params)
    
        
    def plot_fit(self, ax=None, Bode=False):
        
        circuit_list, circuit_funcs = circuits.list_circuits()
        circuitfunc = circuit_funcs[circuit_list.index(self.circuit)]
        
        self.fits = circuitfunc(self.freqs, self.params)
              
        
        if ax is not None:
            if Bode:
                ax.plot(self.freqs, np.abs(self.Z), 'o', color=colors[0])
                ax.plot(self.freqs, np.abs(self.fits), '-', color=colors[0])
                ax2 = ax.twinx()
                ax2.plot(self.freqs, np.angle(self.Z, deg=True), 'x', color=colors[1])
                ax2.plot(self.freqs, np.angle(self.fits, deg=True), '-', color=colors[1])
                ax.set_xscale('log')
                ax.set_xlabel('Frequency/ Hz')
                ax.set_ylabel('|Z|/ $\Omega$')
                ax2.set_ylabel('Phase/ $\degree$')
                
            else:
                ax.plot(self.re/1e6, -self.im/1e6, 'o')
                ax.plot(np.real(self.fits)/1e6, -np.imag(self.fits)/1e6, '-')
                ax.set_xlabel("Z'/ M$\Omega$")
                ax.set_ylabel("Z''/ M$\Omega$")
        




def fit_all_runs(path, circuit, bounds):
    '''
    Fit all specta in a given folder. Assumes spectra are 
    
    named in order, and have only small perturbations between spectra
    
    (i.e. time series of sequential spectra). Takes previous spectrum's
    
    best-fit parameters as initial guess for next spectrum fit.
    
    
    Spectra should be called e.g. 0001s.txt, 
    and have tab-separated format:
    
    <Frequency>	<Re(Z)>	<Im(Z)>
    100.0	19650307.25839891	-11025535.444937838
    110.0	18197021.855608918	-10204987.541515699
    ...
    
    

    Parameters
    ----------
    path : String.
        Path to directory containing multiple EIS spectra.

    Returns
    -------
    d : Dict
        Dictionary of DataFile classes. d[i] is the fit
        of file i.

    '''
    starting_time = time.time()
    
        
    d = {}
    i = 1
    
            
    for f in os.listdir(path):
        
        # Iterate through all files in folder
        if f.endswith('s.txt'):
            
            file = os.path.join(path, f)
                        
            d[i] = DataFile(file, circuit, bounds)
            
            if i == 1:
                # Start new fit routine with genetic algorithm
                n = 0
                fig, ax = plt.subplots()
                while d[i].score > 100000:
                    if n > 5:
                        break
                    d[i].ga_fit(n_iter=100, ax = ax)
                    d[i].LEVM_fit()
                    n += 1
                                    
                print('File 1 fit complete in ', 
                      time.time() - starting_time, 's.')
                

            else:
                # Copy initial parameters from previous fit
                d[i].ga_fit(starting_guess = d[i-1].params,
                            n_iter=10)
                d[i].LEVM_fit()
                    
                if i%50 == 0:
                    print('File %s completed.' %i)


            i += 1
            
                
    # Save params in DataFrame
    df = pd.DataFrame([d[i].params for i in d])
    
    # Save as csv
    out_path = path + '\\fits.csv'
    df.to_csv(out_path)
    
    print('Run complete. Total time: ', 
          time.time() - starting_time, 's.')
    
    return d, df







# if __name__ == '__main__':
    
#     d, df = fit_all_runs(path)
    
    
    
# #%%    Plot parameters vs time
#     param_list = [key for key, item in d[1].params.items()]
    
#     for param in param_list:
#         l = []
#         fig, ax = plt.subplots()
#         for i in d:
#             l.append(d[i].params[param])
#         ax.plot(np.arange(0,45,0.1), l, '.')
#         ax.set_xlabel('Time/ s')
#         ax.set_ylabel(param)
#         ax.set_title(param)
    
    
#     scores = []
#     for i in d:
#         scores.append(d[i].score)
#     fig, ax = plt.subplots()
#     ax.plot(range(len(d)), scores)
#     ax.set_ylabel('Score')



