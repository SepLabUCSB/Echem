import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def hill_equation(C, n, Keq, A, B):
    return A + B * (C**n)/(Keq + C**n)


def backwards_hill_eq(k, n, Keq, A, B):
    a = (k - A)/B
    return (Keq*a/(1-a))**(1/n)


# Vanco
# target = 'Vancomycin'
# file = 'Z:/Projects/Brian/5 - Plaxco collab/20220405 in vivo vanco/meisp_fits_vanco.xlsx'
# A    = 20
# B    = 33.8
# Keq  = 4.16e-4
# n    = 0.85

# # Tobramycin
# target = 'Tobramycin'
# file = 'Z:/Projects/Brian/5 - Plaxco collab/20220406 in vivo tobramycin/meisp_fits_tobramycin.xlsx'
# A    = 102.7
# B    = 966
# Keq  = 4.86e-3
# n    = 0.77

# # Phe
target = 'Phenylalanine'
# file = 'Z:/Projects/Brian/5 - Plaxco collab/20220413 in vivo phenylalanine/meisp_fits_phe.xlsx'
file = 'Z:/Projects/Brian/5 - Plaxco collab/20220427 in vivo phenylalanine fasted/meisp_fits.xlsx'
A    = 62.3
B    = 146
Keq  = 8.66e-2
n    = 0.58



df = pd.read_excel(file)

# df = df[df['ket'] < 140]

# filter_ket = signal.savgol_filter(df['ket'], 51, 1, mode='nearest')

conc = backwards_hill_eq(df['ket'], n, Keq, A, B)
# filt_conc = signal.savgol_filter(conc, 51, 1)
# filt_conc = backwards_hill_eq(filter_ket, n, Keq, A, B)

fig, ax = plt.subplots()
ax.plot(df['time/s']/60, conc)
# ax.plot(df['time/s']/60, filt_conc)
ax.set_xlabel('Time/ min')
ax.set_ylabel(f'[{target}]/ M')
# ax.set_xticks([0,30,60,90,120,150,180,210])