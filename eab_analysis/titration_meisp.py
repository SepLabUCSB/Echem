import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import os

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-07-21'
# data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-08-19 vancomycin bovine blood titration'
# data_dir = r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2022-11-28/Titration'
meisp_dir = data_dir + '-MEISP'

target = 'Vancomycin'

data_dir = r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2023-01-27'
meisp_dir = data_dir + '-MEISP'
target = 'Cocaine'
concs = [
    5e-8,
    1e-7,
    5e-7,
    1e-6,
    5e-5,
    1e-4,
    5e-4,
    1e-3
    ]

co = None # Cut off last co concentrations from all plots



# concs = [0.0E+00,
# 1.0E-07,
# 1.8E-07,
# 3.2E-07,
# 5.6E-07,
# 1.0E-06,
# 1.8E-06,
# 3.2E-06,
# 5.6E-06,
# 1.0E-05,
# 1.8E-05,
# 3.2E-05,
# 5.6E-05,
# 1.0E-04,
# 1.8E-04,
# 3.2E-04,
# 5.6E-04,
# 1.0E-03,
# 1.8E-03,
# 3.2E-03,
# ]

class Spectrum:
        
    def __init__(self, conc, elec, freqs, re, im, params=None,
                 devs = None):
        
        self.conc   = float(conc)
        self.freqs  = freqs.astype(int)
        self.re     = re
        self.im     = im
        self.elec   = elec
        
        self.Z = re + 1j*im
        self.phase = np.angle(self.Z, deg=True)
        self.params = params  # Fit parameters
        self.devs = devs      # Standard devation of fit params
        self.ket = None



def hill_equation(C, n, Keq, A, B):
    # n = 1
    return A + B * (C**n)/(Keq**n + C**n)



def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
            


def extract_data(path):
    
    d = {}
        
    for folder in os.listdir(path):
        
        num, conc = folder.split('_')
        
        conc = concs[int(conc)-1]
        
        f = []
        re = []
        im = []
        
        for file in os.listdir(os.path.join(path, folder)):
            file = os.path.join(path, folder, file)
            if (file.endswith('fits.txt') or
                file.endswith('currents.txt')):
                continue
            if file.endswith('s.txt'):
                fs, res, ims = np.loadtxt(file, unpack=True, skiprows=1)
                f.append(fs)
                re.append(res)
                im.append(ims)
        
        f = np.mean(f, axis=0)
        re = np.mean(re, axis=0)
        im = np.mean(im, axis=0)
        
        if num not in d:
            d[num] = []
            d[num].append(Spectrum(conc, num, f, re, im))
        
        else:
            d[num].append(Spectrum(conc, num, f, re, im))
    
    
    for num in d:
        d[num].sort(key=lambda s:s.conc)        
    
    return d



def write_meisp_files(path, num, specs):
    # Write averaged Z data to MEISP format .txt file
    #specs = list of spectra for an electrode
    
    save_path = os.path.join(path, str(num))
    
    if not os.path.isdir(save_path):
        createFolder(save_path)
    
    for spec in specs:
    
        d = pd.DataFrame(
                {'f': spec.freqs,
                're': spec.re,
                'im': spec.im}
                )
            
        fname = save_path + f'//{spec.conc:.1e}.txt'
    
        d.to_csv(fname, columns = ['f', 're', 'im'],
                     header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], 
                     sep = '\t', index = False, encoding='ascii')
        
        

def extract_meisp_data(path, num, circuit_elements = None):
    # Read fit data from .par and .dev files
    
    path = os.path.join(path, str(num))
    
    par_file = os.path.join(path, f'{num}.par')
    dev_file = os.path.join(path, f'{num}.dev')
    
    
    if not circuit_elements:
        circuit_elements = ['Rs', 'Rct', 'Cdl', 'Cad']
    
    
    # _, _, Rs, Rct, Cdl, Cad = np.loadtxt(par_file, unpack=True)
    # _, _, _, Rs_dev, Rct_dev, Cdl_dev, Cad_dev = np.loadtxt(dev_file, 
    #                                                      unpack=True)
    
    ps = pd.read_fwf(par_file, names=['ind', '_', *circuit_elements])
    ds = pd.read_fwf(dev_file, names=['ind', '_', 's', *circuit_elements])
    
    
    concs = []
    specs = []
    
    
    i = 0
    for file in os.listdir(path):
        if file.endswith('.txt') and len(file) == 11:
            
            concs.append(file)
            
            fs, res, ims = np.loadtxt(os.path.join(path,file), 
                                      unpack=True, skiprows=1)
            
            
            params = {elem: ps[elem][i] for elem in circuit_elements}
            devs   = {elem: ds[elem][i] for elem in circuit_elements}
            
            
            specs.append(Spectrum(float(file.strip('.txt')), 
                                  num, fs, res, ims,
                                  params=params,
                                  devs=devs
                                  )
                          )
         
            i += 1
    
    specs.sort(key=lambda s:s.conc)
    
    for spec in specs:
        spec.ket = 1/(2*spec.params['Rct']*spec.params['Cad'])
        spec.params['ket'] = spec.ket
    
    return specs
 


def plot_phase_maps(specs, name):
    specs = specs[:co]
    colors = plt.cm.magma(np.linspace(0.0,0.7, len(specs)))
    
    fig, ax = plt.subplots()
    i = 0
    for spec in specs:
        ax.plot(spec.freqs, spec.phase, '.-', color = colors[i])
        i += 1
    ax.set_xscale('log')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')
    ax.set_xticks([1e0, 1e1, 1e2, 1e3, 1e4])
    ax.set_xlim(1,1000)
    ax.set_ylim(-85, -65)
    ax.set_yticks([-85, -80, -75, -70, -65])
    # ax.set_title(name)
    
     
    
    
def plot_Z(specs, name):
    specs = specs[:co]
    colors = plt.cm.magma(np.linspace(0,0.7, len(specs)))
    
    fig, ax = plt.subplots()
    i = 0
    for spec in specs:
        ax.plot(spec.freqs, np.absolute(spec.Z)/1000, color = colors[i])
        i += 1
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('|Z|/ k$\Omega$')
    ax.set_xticks([1e0, 1e1, 1e2, 1e3, 1e4])
    ax.set_xlim(1,1000)
    # ax.set_title(name)
            
    

def plot_params(d, title=None):
    
    # plot normalized params vs concentration
    fig, ax = plt.subplots()
    
    names = {'Rs': '$R_{s}$', 
             'Rct': '$R_{ct}$', 
             'Cdl': '$C_{dl}$', 
             'Cad': '$C_{ad}$'}
    
    i = 0
    for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
    # for param in ['Rs', 'Rct']:
        l = []
        concs = []
        for num in d:
            vals = [spec.params[param] for spec in d[num]]
            vals = vals/vals[0]
            l.append(vals)
            concs = [spec.conc for spec in d[num]][:co]
        avg = np.mean(l, axis=0)[:co]
        std = np.std(l, axis=0)[:co]
                
        # from scipy import interpolate
        # x = concs
        # y = avg
        
        # n_interior_knots = 4
        # qs = np.linspace(0, 1, n_interior_knots+2)[1:-1]
        # knots = np.quantile(x, qs)
        # tck = interpolate.splrep(x, y, t=knots, k=3)
        # ys_smooth = interpolate.splev(x, tck)
        
        # ax.plot(x, ys_smooth, 'k-', alpha=0.5, color=colors[i])
    
        ax.errorbar(concs, avg, std, marker='o', linestyle='',
                    capsize=4, elinewidth=2,
                    label= names[param], color=colors[i])
        # ax.plot(xnew, f_cubic(xnew), color=colors[i])
        i += 1
        
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('Normalized Parameter')
    ax.set_xscale('log') 
    ax.set_xticks([1e-7,1e-6,1e-5,1e-4,1e-3])
    ax.set_yticks([0.4,0.6,0.8,1.0])
    ax.set_ylim(ax.get_ylim()[0], 1.1)
    if title:
        ax.set_title(title)
    plt.legend()
        
    
    
    
    # plot ket vs concentration
    kets = []
    concs = []
    for num in d:
        ket = np.array([spec.ket for spec in d[num]])
        kets.append(ket)
        concs = [spec.conc for spec in d[num]][:co]
    avg = np.mean(kets, axis=0)[:co]
    std = np.std(kets, axis=0)[:co]
    
    fig, ax = plt.subplots()
    ax.errorbar(concs, avg, std, capsize=4, 
                    elinewidth=2,)
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.set_xscale('log') 
    ax.set_xticks([1e-7,1e-6,1e-5,1e-4,1e-3])
    if title:
        ax.set_title(title)
    


def save_fits(d, data_dir):  
    
    ls = []
    for num in d:
        for spec in d[num]:
            
            l = [spec.elec,
                 spec.conc,
                 *[val for elem, val in spec.params.items()],
                 *[dev for elem, dev in spec.devs.items()]
                 ]
            ls.append(l)
    
    names = ['elec',
             'conc',
             *[elem for elem in spec.params],
             *[f'{elem}_dev' for elem in spec.devs]]
    
    df = pd.DataFrame(ls, columns = names)
    
    writer = pd.ExcelWriter(data_dir + '/fits.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {data_dir}/fits.xlsx')



def hill_fit(fit_d):
    kets = []
    concs = []
    for num in fit_d:
        ket = np.array([spec.ket for spec in fit_d[num]])
        kets.append(ket)
        concs = np.array([spec.conc for spec in fit_d[num]])[:co]
    
    
    avg = np.mean(kets, axis=0)[:co]
    std = np.std(kets, axis=0)[:co]
        
    
    popt, pcov = optimize.curve_fit(hill_equation, concs, avg, sigma=std,
                                          maxfev=100000)
    
    
    s1 = [r'$k_{et}$ = $k_{min}$ + $\Delta$k '+ r'$\frac{C^{n}}{K_{D}^{n} + C^{n}}$',
          '',
          r'$k_{min}$  = ' + f'{popt[2]:.1f}' + ' $s^{-1}$',
          r'$\Delta$k    = ' + f'{popt[3]:.1f}' + ' $s^{-1}$',
          r'$K_{D}$' + f'    = {popt[1]/1e-6:.0f} $\mu$M',
          f'n      = {popt[0]:.2}',
    ]
    
    s = '\n'.join([ss for ss in s1])
    
    fig, ax = plt.subplots()
    ax.errorbar(concs, avg, std, capsize=4, elinewidth=2, ecolor='k',
                marker='o', color='k', lw=0)
    ax.plot(concs, hill_equation(concs, *popt), 'red')
    ax.set_xscale('log')
    ax.set_xticks([1e-7,1e-6,1e-5,1e-4,1e-3])
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.text(2e-7, 55, s)
    
    return popt, pcov


#  Run functions        
   
circuit_elements = ['Rs', 'Rct', 'Cad', 'phi', 'Cdl']
     
d = extract_data(data_dir)
fit_d = {}


for num in d:
    write_meisp_files(meisp_dir, num, d[num])
    fit_d[num] = extract_meisp_data(meisp_dir, num, circuit_elements)
    plot_phase_maps(d[num], f'Electrode {num}')
    # plot_Z(d[num], f'Electrode {num}')


plot_params(fit_d)
# save_fits(fit_d, meisp_dir)
popt, pcov = hill_fit(fit_d)



















