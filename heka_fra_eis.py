import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

from io import StringIO

plt.style.use('Z:\Projects\Brian\scientific.mplstyle')







def extract_data(file):
    s = StringIO()
        
    def isfloat(x):
        try:
            float(x)
            return True
        except:
            return False
    
    with open(file, 'r') as f:
        for line in f:
            if isfloat(line.split(',')[0]):
                #skip rows which don't start with the index number (1, 2, ...)
                s.write(line)
    
    s.seek(0)
    df = pd.read_csv(s, names=(
            'f','logf','AY', 'logAY','RY','IY','phase','RZ','IZ','modZ', 'logmodZ'
            ), engine='c', dtype=float)
    
    f = df['f'].to_numpy()
    f = np.around(f, 3)
    re = df['RZ'].to_numpy()
    im = -df['IZ'].to_numpy()    
    
    return f, re, im


def save_corr_file(file, f, Z_corr, phase_corr):
    folder = 'C:/Users/BRoehrich/Desktop/git/echem/fra_eis_corrections'
    path = os.path.join(folder, file)
        
    df = pd.DataFrame({'f': f, 'Z_corr':Z_corr, 'phase_corr':phase_corr})
    df.to_csv(path, index=False)
    return


def get_corrections(R, f, re, im):
    Z = np.absolute(re + 1j*im)
    phase = np.angle(re + 1j*im, deg=True)
    
    Z_corr = Z/R
    phase_corr = -phase
    
    return f, Z_corr, phase_corr


def make_correction_file(R, f, re, im, file):
    f, Z_corr, phase_corr = get_corrections(R,f,re,im)
    save_corr_file(file,f,Z_corr,phase_corr)
    return


def correct_Z(f, re, im, corr_file='C:/Users/BRoehrich/Desktop/git/echem/fra_eis_corrections/10M.csv'):
    corr_df = pd.read_csv(corr_file)
    corr_df['f'] = np.around(corr_df['f'].to_numpy(), 3)
           
    Z = np.absolute(re + 1j*im)
    phase = np.angle(re + 1j*im, deg=True)

    df = pd.DataFrame({'f':f, 'Z': Z, 'phase':phase})
    if all(df['f'] == corr_df['f']):
        df['Z'] /= corr_df['Z_corr']
        df['phase'] += corr_df['phase_corr']
    else:
        print(f'Frequencies dont match correction file {corr_file}')
        print(list(zip(df['f'], df['f'] == corr_df['f'])))
        print('Not correcting')
            
    
    Z = df['Z'] * np.exp(1j*df['phase']*np.pi/180)
    
    return Z
    
        
def plot_Nyquist(Z, dpi=300):
    fig, ax = plt.subplots(dpi=dpi)
    ax.plot(np.real(Z)/1e6, -np.imag(Z)/1e6, 'o-')
    ax.set_xlabel(r"Z' / M$\Omega$")
    ax.set_ylabel(r"Z'' / M$\Omega$")
    
    low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    

def plot_Bodes(f, Z, dpi=300):
    
    phase = np.angle(Z, deg=True)
    mod   = np.absolute(Z)
    
    fig, ax = plt.subplots(dpi=dpi)
    ax.plot(f, phase, 'o-')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')
    ax.set_xscale('log')
    
    fig, ax = plt.subplots(dpi=dpi)
    ax.plot(f, mod/1e6, 'o-')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('|Z|/ M$\Omega$')
    ax.set_xscale('log')
    
def save_Z(f, Z, fout):
    d = pd.DataFrame(
            {'f': f,
            're': np.real(Z),
            'im': np.imag(Z)}
            )
        
    
    if not fout.endswith('.txt'):
        fout = fout.split('.')[0]
        fout += '.txt'
            
    
    d.to_csv(fout, columns = ['f', 're', 'im'],
                 header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], 
                 sep = '\t', index = False, encoding='ascii')
    print(f'Saved as {fout}')
    
    

# files = ["D:/Brian/nanoelectrode tests/20221028/Pt nanoE 1mM FcMeOH 100mM KCl 10k to 100m.asc"]  
# files = ['D:/Brian/nanoelectrode tests/20221028/Pt nanoE 1mM K4FeCN6 100mM KCl 10k to 100m.asc']  
files = ['D:/Brian/nanoelectrode tests/20221031/1mM FcMeOH 05 mM KCl Pt UME 10-19-10.asc',
         'D:/Brian/nanoelectrode tests/20221031/1mM FcMeOH 05 mM KCl Pt UME 10-19-7.asc']
# files = ['D:/Brian/nanoelectrode tests/20221027 eis tests/10M 10k to 100m_for calibration.asc']

# fig, ax = plt.subplots()

for file in files:
    f, re, im = extract_data(file)
    Z = correct_Z(f, re, im)
    plot_Nyquist(Z)
    save_Z(f, Z, file)

# ax.set_xscale('log')



   