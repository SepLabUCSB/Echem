import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import os

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

data_dir = r'C:\Users\BRoehrich\Desktop\2022-05-10\2022-05-10'
meisp_dir = r'C:\Users\BRoehrich\Desktop\2022-05-10\2022-05-10-MEISP'

target = 'Phenylalanine'

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
    return A + B * (C**n)/(Keq + C**n)



def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
            


def extract_data(path):
    
    d = {}
        
    for folder in os.listdir(path):
        
        num, mv, conc = folder.split('_')
        num  = num + '_' + mv
        conc = float(conc)
        
        f = []
        re = []
        im = []
        
        for file in os.listdir(os.path.join(path, folder)):
            file = os.path.join(path, folder, file)
            if file.endswith('fits.txt'):
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
            
            # specs.append(Spectrum(float(file.strip('.txt')), 
            #                       num, fs, res, ims,
            #                       params={'Rs': Rs[i],
            #                               'Rct': Rct[i],
            #                               'Cdl': Cdl[i],
            #                               'Cad': Cad[i]},
            #                       devs={'Rs': Rs_dev[i],
            #                               'Rct': Rct_dev[i],
            #                               'Cdl': Cdl_dev[i],
            #                               'Cad': Cad_dev[i]}
            #                       )
            #              )
        
            i += 1
    
    specs.sort(key=lambda s:s.conc)
    
    for spec in specs:
        spec.ket = 1/(2*spec.params['Rct']*spec.params['Cad'])
        spec.params['ket'] = spec.ket
    
    return specs
 


def plot_phase_maps(specs, name):
    
    colors = plt.cm.viridis(np.linspace(0.1,0.8, len(specs)))
    
    fig, ax = plt.subplots()
    i = 0
    for spec in specs:
        ax.plot(spec.freqs, spec.phase, color = colors[i])
        i += 1
    ax.set_xscale('log')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')
    ax.set_xticks([1e0, 1e1, 1e2, 1e3, 1e4])
    ax.set_title(name)
            

def plot_params(d, title=None):
    
    # plot normalized params vs concentration
    fig, ax = plt.subplots()
    
    names = {'Rs': '$R_{s}$', 
             'Rct': '$R_{ct}$', 
             'Cdl': '$C_{dl}$', 
             'Cad': '$C_{ad}$'}
    
    for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
        l = []
        concs = []
        for num in d:
            vals = [spec.params[param] for spec in d[num]]
            vals = vals/vals[0]
            l.append(vals)
            concs = [spec.conc for spec in d[num]]
        avg = np.mean(l, axis=0)
        std = np.std(l, axis=0)
    
        ax.errorbar(concs, avg, std, capsize=4, elinewidth=2, 
                    label= names[param])
    
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('Normalized Parameter')
    ax.set_xscale('log') 
    if title:
        ax.set_title(title)
    plt.legend()
        
    
    
    
    # plot ket vs concentration
    kets = []
    concs = []
    for num in d:
        ket = np.array([spec.ket for spec in d[num]])
        kets.append(ket)
        concs = [spec.conc for spec in d[num]]
    avg = np.mean(kets, axis=0)
    std = np.std(kets, axis=0)
    
    fig, ax = plt.subplots()
    ax.errorbar(concs, avg, std, capsize=4, 
                    elinewidth=2,)
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.set_xscale('log') 
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
    print(df)
    
    
    
    
    # df = pd.DataFrame(array, columns=names)
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
        concs = np.array([spec.conc for spec in fit_d[num]])
    
    
    avg = np.mean(kets, axis=0)
    std = np.std(kets, axis=0)
    
    # avg, std, concs = (avg[:-3], std[:-3], concs[:-3])
    
    
    popt, pcov = optimize.curve_fit(hill_equation, concs, avg, sigma=std,
                                          maxfev=10000)
    
    
    s1 = [r'$k_{et}$ = A + B '+ r'$\frac{C^{n}}{K_{D} + C^{n}}$',
          '',
          f'A    = {popt[2]:.1f}' + ' $s^{-1}$',
          f'B    = {popt[3]:.1f}' + ' M $s^{-1}$',
          r'$K_{D}$' + f'  = {popt[1]:.2e} M',
          f'n    = {popt[0]:.2}',
    ]
    
    s = '\n'.join([ss for ss in s1])
    
    fig, ax = plt.subplots()
    ax.errorbar(concs, avg, std, capsize=4, elinewidth=2, ecolor='k',
                marker='o', color='k', lw=0)
    ax.plot(concs, hill_equation(concs, *popt), 'red')
    ax.set_xscale('log')
    ax.set_xlabel(f'[{target}]/ M')
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.text(2e-7, 90, s)
    
    return popt, pcov


#%%  Run functions        
   
circuit_elements = ['Rs', 'Rct', 'Cad', 'phi', 'Cdl']
     
d = extract_data(data_dir)
# d = [1,2,3,4]
fit_d = {}


for num in d:
    # write_meisp_files(meisp_dir, num, d[num])
    fit_d[num] = extract_meisp_data(meisp_dir, num, circuit_elements)
    # plot_phase_maps(d[num], f'Electrode {num}')


# for num in fit_d:
#     fit_d[num] = fit_d[num][:-2]

# plot_params(fit_d)
save_fits(fit_d, meisp_dir)
# popt, pcov = hill_fit(fit_d)



















