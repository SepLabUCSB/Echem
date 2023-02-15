import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal, optimize
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


## Vanco
# 4/5
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20220405 in vivo vanco/meisp_fits_vanco.xlsx'
# 7/27
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20220727 vanco invivo/meisp_fits.xlsx'
# 9/14
file = r'Z:/Projects/Brian/5 - Plaxco collab/20220914 vanco invivo 25mV good/meisp_fits.xlsx'
# 10/6 rat 1
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20221005 vanco invivo 25mV 2 rats/rat 1 meisp_fits.xlsx'
# 10/6 rat 2
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20221005 vanco invivo 25mV 2 rats/rat 2 meisp_fits.xlsx'

# 11/18 fast in-vitro recording
# file = r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2022-11-18 fast recording vanco invitro/meisp_fits.xlsx'
# 11/28 fast PBS recording
# file = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-11-28\fast recording additions/meisp_fits.xlsx'


## Phe
# 4/13
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20220413 in vivo phenylalanine/meisp_fits_phe.xlsx'
# 4/27
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20220427 in vivo phenylalanine fasted/meisp_fits.xlsx'
# 9/20
# file = r'Z:/Projects/Brian/5 - Plaxco collab/20220920 phe invivo 25mVAC/meisp_fits.xlsx'



# def f(k, n, Keq, A, B):
#     a = (k-A)/B
#     return ((a*Keq)/(1-a))**(1.0/n)


# def f(k, n, KD, A, B):
#     return ((KD*(k-A))/(A+B-k))**(1/n)

def f(k, n, KD, A, B):
    a = (k - A)/B
    return ((a*KD**n)/(1-a))**(1/n)


def hill_equation(C, n, Keq, A, B):
    # n = 1
    return A + B * ((C**n)/(Keq**n + C**n))


def fit_decay(t, conc, add=0):
        
    smoothed = signal.savgol_filter(conc, 101, 1, mode='constant')
    peak = np.where(smoothed == max(smoothed))[0][0]
    # fig, ax = plt.subplots()
    # ax.plot(conc)
    # ax.plot(smoothed)
    
    peak += add
    
    ts = t[peak:] - t[peak]
    ys = conc[peak:]
    
    
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        
    popt, pcov = optimize.curve_fit(decay_eq, ts, ys, maxfev=100000)
    
    plot_ts = ts + t[peak]
    plot_ys = decay_eq(ts, *popt)
    
    # fig, ax = plt.subplots()
    # ax.plot(t, conc)
    # ax.plot(plot_ts, plot_ys)
    
    return plot_ts, plot_ys, popt, pcov



def plot_params(df):
    fig, ax = plt.subplots()
    for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
        ax.plot(df['t']/60, df[param]/df[param][0], label=param)
    
    ax.set_xlabel('Time/ min')
    ax.set_ylabel('Normalized parameter')
    ax.set_xticks(np.arange(0, ax.get_xlim()[1], 60))
    # ax.set_xlim(ax.get_xlim()[0], 250)
    # ymin, ymax = ax.get_ylim()
    # ax.set_xlim(ax.get_xlim()[0],
    #             1.2*ax.get_xlim()[1])
    # ax.set_ylim(0.3, ymax)
    ax.legend()
    


# Vancomycin
n, KD, A, B = [1., 1.44182502e-04, 2.39721831e+01, 4.31584165e+01] # Normal
# n, KD, A, B = [1.17611897e+00, 8.60122501e-05, 2.62143735e+01, 3.94513993e+01] # cutoff 10Hz to 1kHz fit
# n, KD, A, B = [0.77, 197e-6, 31.1, 66.3] # True 10Hz to 1kHz

# Phe...
# n, KD, A, B = [3.78042233e-01, 6.89266567e-03, 6.29065917e+01, 1.62870307e+02]


# t,Rs,Rct,Cdl,Cad,ket,_,_,_,_ = np.loadtxt(file, unpack=True, skiprows=1)
df = pd.read_excel(file, skiprows=1, names=('t','Rs','Rct','Cad','phi','Cdl', 'ket',
                                            'a', 'b', 'c', 'd', 'e'))

# df = pd.read_excel(file, skiprows=1, names=('t','Rs','Rct','Cdl','Cad', 'ket',
#                                             'a', 'b', 'c', 'd'))


df = df.dropna()

def set_values(value, arr, ref_arr, lim1, lim2):
    idxs = [i for i, val in enumerate(ref_arr) if (val >= lim1 and val < lim2)]
    for i, _ in enumerate(arr):
        if i in idxs:
            arr[i] = value
    return arr


t = df['t'].to_numpy()
Rs = df['Rs'].to_numpy()
Rct = df['Rct'].to_numpy()
Cdl = df['Cdl'].to_numpy()
Cad = df['Cad'].to_numpy()
ket = df['ket'].to_numpy()

conc = f(ket, n, KD, A, B)
# conc = 10**((ket-184.1)/23.3)

#%%

# arr = np.ones(len(conc))

# arr = set_values(0, arr, t/60, 0, 5)
# arr = set_values(1e-6, arr, t/60, 5, 11)
# arr = set_values(4.5e-5, arr, t/60, 11, 16)
# arr = set_values(9.5e-5, arr, t/60, 15.1, 1000)
# applied_conc = arr

fig, ax = plt.subplots(figsize=(5,5), dpi=300)
ax.plot(t/60, conc/1e-6, '.', color = colors[3], alpha=0.1)
# ax.plot(t/60, applied_conc, '--', color='k')
# ax.plot(t/60, signal.savgol_filter(conc, 7, 1), color=colors[3])
ax.set_xlabel('Time/ min')
ax.set_ylabel(r'[Vancomycin]/ M')
# ax.set_yscale('log')
ax.set_ylim(-10, 105)
# ax.text(12, 4e-6,'Applied conc', color='k')
# ax.text(10, 7e-4, 'Measured conc', color=colors[2])

#%%
# plot_params(df)


# def decay_eq(t, t0, A, alpha, C):
#     return A*np.exp(-(1/alpha)*(t-t0)) + C

def decay_eq(t, t0, A, alpha, B, beta, C):
    return A*np.exp(-(1/alpha)*(t-t0)) + B*np.exp(-(1/beta)*(t-t0)) + C




# half = int(len(t)/2)
# plot_ts, plot_ys, popt, pcov = fit_decay(t[:half], conc[:half],1)
# # plot_ts, plot_ys, popt, pcov = fit_decay(t, conc,-5)

# # t0, A, alpha, C = popt
# # _, _, std, _ = np.sqrt(np.diag(pcov))

# t0, A, alpha, B, beta, C = popt
# _, _, std, _, std1, _ = np.sqrt(np.diag(pcov))

# # alpha = 0
# # std = 0


# print(f'Decay constant: {alpha/60:0.2f} min')
# s1 = r'$\tau_{1}$ : '
# s2 = f'{alpha/60:0.2f} '
# s3 = r'$\pm$'
# s4 = f' {std/60:0.2f}'
# s = s1+s2+s3+s4

# s1 = r'$\tau_{2}$ : '
# s2 = f'{beta/60:0.2f} '
# s3 = r'$\pm$'
# s4 = f' {std1/60:0.2f}'
# s2 = s1+s2+s3+s4

# s = s + '\n' + s2




# fig, ax = plt.subplots(figsize=(5,5), dpi=300)
# ax.plot(t/60, conc/1e-6, '.', color = colors[2], alpha=0.1)
# ax.plot(t/60, signal.savgol_filter(conc/1e-6, 7, 1), color=colors[2])
# ax.plot(plot_ts/60, plot_ys/1e-6, color= 'k')
# ax.text(72,400, s)
# # ax.set_ylim(-11, 102)



# plot_ts, plot_ys, popt, pcov = fit_decay(t[half:], conc[half:], add=0)


# # t0, A, alpha, C = popt
# # _, _, std, _ = np.sqrt(np.diag(pcov))

# t0, A, alpha, B, beta, C = popt
# _, _, std, _, std1, _ = np.sqrt(np.diag(pcov))



# print(f'Decay constant: {alpha/60:0.2f} min')
# s1 = r'$\tau_{1}$ : '
# s2 = f'{alpha/60:0.2f} '
# s3 = r'$\pm$'
# s4 = f' {std/60:0.2f}'
# s = s1+s2+s3+s4

# s1 = r'$\tau_{2}$ : '
# s2 = f'{beta/60:0.2f} '
# s3 = r'$\pm$'
# s4 = f' {std1/60:0.2f}'
# s2 = s1+s2+s3+s4

# s = s + '\n' + s2


# ax.plot(plot_ts/60, plot_ys/1e-6, color= 'k')
# ax.set_xlabel('Time/ min')
# ax.set_ylabel('[Phenylalanine]/ $\mu$M')
# ax.set_xticks(np.arange(0, ax.get_xlim()[1], 60))
# ax.text(72, 300, s)

# ax.set_ylim(-50, 1000)


# fig, ax = plt.subplots(figsize=(5,5), dpi=300)
# ax.plot(t/60, ket, '.', color = colors[0])
# ax.set_xlabel('Time/ min')
# ax.set_ylabel('$k_{et}$/ $s^{-1}$')
# ax.set_xticks(np.arange(0, ax.get_xlim()[1], 60))