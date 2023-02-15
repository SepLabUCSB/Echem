import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
phe_fit_file = r"Z:/Projects/Brian/5 - Plaxco collab/20220615 Phe rat's blood titration/no CPE fits.xlsx"
invivo_fits = r'Z:/Projects/Brian/5 - Plaxco collab/20220920 phe invivo 25mVAC/meisp_fits.xlsx'



def hill_equation(C, n, Keq, A, B):
    # n = 1
    return A + B * (C**n)/(Keq**n + C**n)


df = pd.read_excel(phe_fit_file)

concs = df[df['elec'] == 1]['conc'].to_numpy()
ket = np.average([ df[df['elec'] == i]['ket'].to_numpy() for i in [1,2,3,4]], axis=0)
ket_std = np.std([ df[df['elec'] == i]['ket'].to_numpy() for i in [1,2,3,4]], axis=0)


x = concs[20:-5]
y = ket[20:-5]
ystd = ket_std[20:-5]


bounds = (
          np.array([0, 1e-5, -np.inf, -np.inf]),
          np.array([np.inf, 100, np.inf, np.inf])
          )

popt, pcov = optimize.curve_fit(hill_equation, x, y, sigma=ystd, maxfev=100000,
                                bounds=bounds)

s1 = [r'$k_{et}$ = $k_{min}$ + $\Delta$k '+ r'$\frac{C^{n}}{K_{D}^{n} + C^{n}}$',
          '',
          r'$k_{min}$  = ' + f'{popt[2]:.1f}' + ' $s^{-1}$',
          r'$\Delta$k    = ' + f'{popt[3]:.1f}' + ' $s^{-1}$',
          r'$K_{D}$' + f'    = {popt[1]/1e-6:.0f} $\mu$M',
          f'n      = {popt[0]:.2}',
    ]   
s = '\n'.join([ss for ss in s1])

fig, ax = plt.subplots(dpi=600)
ax.errorbar(x, y, yerr=ystd, capsize=4, elinewidth=2, ecolor='k',
                marker='o', color='k', lw=0)
# ax.errorbar(concs, ket, yerr=ket_std, capsize=4, elinewidth=2, ecolor='k',
#                 marker='o', color='k', lw=0)

# wc = np.logspace(-9, 4, 1000)
ax.plot(x, hill_equation(x, *popt), color='red')
ax.set_xscale('log')
ax.set_xlabel('[Phenylalanine]/ M')
ax.set_ylabel('$k_{et}$/ $s^{-1}$')
# ax.text(4e-5, 100, s)

fig, ax = plt.subplots()
ax.text(0.1, 0.1, s)


