import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from EIS_Fit import EIS_fit

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


data_dir = r'C:/Users/BRoehrich/Desktop/2022-02-23_vancomycin'

target = 'Vancomycin'

starting_bounds = {
    'R1': [1e-1, 1e9],
    'R2': [2000, 300000],
    'Q1': [1e-15, 1],
    'n1': [0.8,1.1],
    'Q2': [1e-15, 1],
    'n2': [0.8,1.2]
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
    'R2': 100000, 
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


class Spectrum:
        
    def __init__(self, C, freqs, re, im, file, AC = None,
                 DC = None, elec_size = None, params = None):
        
        self.file   = file
        self.C      = float(C)
        self.freqs  = freqs.astype(int)
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
        

# def extract_data(folder, d):
#     elec, C = folder.split('_')
    
#     freqs = np.array([])
#     re_list = []
#     im_list = []
    
#     for file in os.listdir(os.path.join(data_dir, folder)): 
#         file = os.path.join(data_dir, folder, file)
#         if file.endswith('s.txt'):
#             df = pd.read_csv(file, skiprows=1, names=('f', 're', 'im'), sep='\t')

#             freqs = df['f'].to_numpy()
#             re_list.append(df['re'].to_numpy())
#             im_list.append(df['im'].to_numpy())
            
            
        
#     re = np.mean(re_list, axis=0)
#     im = np.mean(im_list, axis=0)
    
#     d.append(Spectrum(C, freqs, re, im, file))
    
#     return d



def fit_all(d):

    for i in range(len(d)):
        if i == 0:
            d[i].fit(n_iter = 200, plot=True, starting_guess = starting_guess)
        elif i > 0:
            d[i].fit(n_iter = 50, 
                     starting_guess = {
                         param:np.random.uniform(0.98,1.02)*val 
                         for param, val in d[i-1].params.items()
                                       },
                     bounds = starting_bounds,
                     
                      # bounds = {
                      #     param:[val/3, val*3] for 
                      #     param, val in d[i-1].params.items()
                      #     }, 
                      plot = True)


l = []

electrode_numbers = set()

for folder in os.listdir(data_dir):
    elec, conc = folder.split('_')
    electrode_numbers.add(elec)
    conc = float(conc)
    
    freqs = np.array([])
    re_list = []
    im_list = []  
    
    for file in os.listdir(os.path.join(data_dir, folder)): 
        file = os.path.join(data_dir, folder, file)
        if file.endswith('s.txt'):
            if file.endswith('fits.txt'): continue
            df = pd.read_csv(file, skiprows=1, names=('f', 're', 'im'), sep='\t')

            freqs = df['f'].to_numpy()[:-6]
            re_list.append(df['re'].to_numpy()[:-6])
            im_list.append(df['im'].to_numpy()[:-6])
            
        
    re = np.mean(re_list, axis=0)
    im = np.mean(im_list, axis=0)
    
    spec = Spectrum(conc, freqs, re, im, file, elec_size=elec)
    # spec.fit(n_iter = 100, plot=True, starting_guess = starting_guess)
    if conc < 4e-3:
        l.append(spec)
        
#%% 

def plot_Bodes(l, name=None):
    # Bode plots  
    # l: list of spectra sorted by target concentration
        
    colors = plt.cm.viridis(np.linspace(0.2,0.8, len(l)))
    
    i = 0
    fig, ax = plt.subplots()
    for spectrum in l:
        ax.plot(spectrum.freqs, spectrum.phase, color = colors[i])
        i += 1
    
    ax.set_xscale('log')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')
    if name:
        ax.set_title(name)


def plot_phase_change(l):
    # Plot % phase change vs target concentration
    # l: list of spectra sorted by target concentration
    
    conc     = [s.C for s in l]
    phase4  = [s.phase[1] for s in l]
    phase64 = [s.phase[7] for s in l]
    
    p4_0 = phase4[0]
    p64_0 = phase64[0]
    
    for i in range(len(phase4)):
        phase4[i] = 100*(1-phase4[i]/p4_0)
        phase64[i] = 100*(1-phase64[i]/p64_0)
    
    
    fig, ax = plt.subplots()
    ax.plot(conc, phase4, label='4 Hz')
    ax.plot(conc, phase64, label='64 Hz')
    ax.set_xscale('log')
    ax.set_xlabel('Vancomycin/ M')
    ax.set_ylabel('% Phase Change')
    ax.legend()
    # ax.set_xticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3])
    # ax.set_xlim(1e-7, 1e-3)



def plot_fit_params(d, name):
    # Plot fit parameters vs target concentration
    
    R1s  = np.array([spectrum.params['R1'] for spectrum in d])
    R2s  = np.array([spectrum.params['R2'] for spectrum in d])
    Cdls = np.array([spectrum.params['Q1'] for spectrum in d])
    ndls = np.array([spectrum.params['n1'] for spectrum in d])
    Cas  = np.array([spectrum.params['Q2'] for spectrum in d])
    nas  = np.array([spectrum.params['n2'] for spectrum in d])
    ket  = 1/(2*R2s*Cas)
    
    
    concs = np.array([spectrum.C for spectrum in d])
    
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    # Rs
    # fig, ax = plt.subplots()
    # ax.plot(concs, R1s, color = colors[0])
    # ax.set_xlabel(f'{target}/ M')
    # ax.set_ylabel('$R_{s}$ / $\Omega$')
    # ax.set_xscale('log')
    # ax.set_title(name)
    
    
    # Rct
    fig, ax = plt.subplots()
    ax.plot(concs, R2s/1000, color = colors[1])
    ax.set_xlabel(f'{target}/ M')
    ax.set_ylabel('$R_{ct}$ / $k\Omega$')
    ax.set_xscale('log')
    ax.set_title(name)
    
    
    # Cdl
    fig, ax = plt.subplots()
    ax.plot(concs, Cdls/1e-9, color = colors[2])
    ax.set_xlabel(f'{target}/ M')
    ax.set_ylabel('$C_{dl}$ / nF')
    ax.set_xscale('log')
    ax.set_title(name)
    
    # Ca
    fig, ax = plt.subplots()
    ax.plot(concs, Cas/1e-6, color = colors[3])
    ax.set_xlabel(f'{target}/ M')
    ax.set_ylabel('$C_{ad}$ / $\mu$F')
    ax.set_xscale('log')
    ax.set_title(name)

    # # All 4 parameters normalized
    fig, ax = plt.subplots()
    # ax.plot(concs, R1s/R1s[1], color = colors[0], label='$R_{s}$')
    ax.plot(concs, R2s/R2s[1], color = colors[1], label='$R_{ct}$')
    ax.plot(concs, Cdls/Cdls[1], color = colors[2], label='$C_{dl}$')
    ax.plot(concs, Cas/Cas[1], color = colors[3], label='$C_{ad}$')
    # ax.plot(concs, ndls/ndls[1], color=colors[4], label='$n_{dl}$')
    # ax.plot(concs, nas/nas[1], color=colors[5], label='$n_{ad}$')
    ax.set_xlabel(f'{target}/ M')
    ax.set_ylabel('Normalized Parameter')
    ax.set_xscale('log')
    ax.legend()
    ax.set_title(name)
    
    # ket
    fig, ax = plt.subplots()
    ax.plot(concs, ket, color = colors[4])
    ax.set_xlabel(f'{target}/ M')
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.set_xscale('log')
    ax.set_title(name)
    
    l = [concs, ket, R1s, R2s, Cdls, Cas, ndls, nas]
    
    return l


ds = [[d for d in l if d.elec_size == str(i)] for i in electrode_numbers]


# d1 = [d for d in l if d.elec_size == '1']
# d2 = [d for d in l if d.elec_size == '2']
# d7 = [d for d in l if d.elec_size == '7']
# d8 = [d for d in l if d.elec_size == '8']

for d in ds:
    d.sort(key=lambda d:d.C)
    fit_all(d)
#%%

ls = [plot_fit_params(d, f'Electrode {d[0].elec_size}') for d in ds]
_ = [plot_Bodes(d, f'Electrode {d[0].elec_size}') for d in ds]

# l1 = plot_fit_params(d1, 'Electrode 1')
# l4 = plot_fit_params(d4, 'Electrode 4')
# l5 = plot_fit_params(d5, 'Electrode 5')
# l6 = plot_fit_params(d6, 'Electrode 6')
#%%
kets = [l[1] for l in ls]

mean = np.mean(kets, axis=0)
std = np.std(kets, axis=0)

fig, ax = plt.subplots()
ax.errorbar(ls[0][0], mean, std, capsize=4, 
                    elinewidth=2,)
ax.set_xlabel(f'{target}/ M')
ax.set_ylabel('$k_{et}$/ $s^{-1}$')
ax.set_xscale('log') 
# mean = np.mean([ket1, ket2, ket3], axis=0)
# std = np.std([ket1, ket2, ket3], axis=0)

#%%

names = ['concs', '$k_{et}$', '$R_{s}$', '$R_{ct}$', '$C_{dl}$', 
          '$C_{ad}$', '$n_{dl}$', '$n_{ad}$']

fig, ax = plt.subplots()

for i in range(len(ls[0])):
    if i > 2 and i < 6:
        mean = np.mean([l[i]/l[i][1] for l in ls], axis=0)
        # mean = np.mean([l1[i]/l1[i][1], 
        #                 l4[i]/l4[i][1], 
        #                 l5[i]/l5[i][1],
        #                 l6[i]/l6[i][1]
        #                 ], axis=0)
        std = np.std([l[i]/l[i][1] for l in ls], axis=0)
        # std = np.std([l1[i]/l1[i][1], 
        #                 l4[i]/l4[i][1], 
        #                 l5[i]/l5[i][1],
        #                 l6[i]/l6[i][1]
        #                 ], axis=0)
        ax.errorbar(ls[0][0], mean, std, capsize=4, 
                    elinewidth=2, label= names[i])


# ax.plot(concs, ket1, 'o-', label = 'Electrode 1')
# ax.plot(concs, ket2, 'o-', label = 'Electrode 2')
# ax.plot(concs, ket3, 'o-', label = 'Electrode 3')
# ax.errorbar(concs, mean, std, capsize=3)
ax.set_xlabel(f'{target}/ M')
ax.set_ylabel('Normalized Parameter')
ax.set_xscale('log')
ax.legend()


# ket_mean = np.mean([l1[1], l2[1], l3[1]], axis=0)

# ket_std = np.std([l1[1], l2[1], l3[1]], axis=0)

# fig, ax = plt.subplots()
# ax.errorbar(l1[0], ket_mean, ket_std, capsize=4, 
#                     elinewidth=2,)
# ax.set_xlabel(f'{target}/ M')
# ax.set_ylabel('$k_{et}$/ $s^{-1}$')
# ax.set_xscale('log')


def save_fits(l):
    '''
    l = list of Spectrum objects
    '''
        
    rows_list = []
    
    for spectrum in l:
        dict1 = {
               'Conc': spectrum.C,
               'elec_size': spectrum.elec_size,
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
        rows_list.sort(key= lambda d: (d['elec_size'], d['Conc']))
        
        # out_df = out_df.append(data_list, ignore_index=True)
        df = pd.DataFrame(rows_list)
    
    df.to_csv(r'C:\Users\BRoehrich\Desktop\saved_fits.csv')
    
    return df
        
        
        

        
        
        



        
                