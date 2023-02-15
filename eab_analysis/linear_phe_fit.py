import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
phe_fit_file = r"Z:/Projects/Brian/5 - Plaxco collab/20220615 Phe rat's blood titration/no CPE fits.xlsx"
invivo_fits = r'Z:/Projects/Brian/5 - Plaxco collab/20220920 phe invivo 25mVAC/meisp_fits.xlsx'



df = pd.read_excel(phe_fit_file)

concs = df[df['elec'] == 1]['conc'].to_numpy()
ket = np.average([ df[df['elec'] == i]['ket'].to_numpy() for i in [1,2,3,4]], axis=0)
ket_std = np.std([ df[df['elec'] == i]['ket'].to_numpy() for i in [1,2,3,4]], axis=0)


x = concs[20:-5]
y = ket[20:-5]
ystd = ket_std[20:-5]


fit, cov = np.polyfit(np.log10(x), y, deg=1, cov=True)

s = f'k = {fit[0]:0.2f}*log([Phe]) + {fit[1]:0.2f}'


fig, ax = plt.subplots(dpi=600)
ax.errorbar(x, y, yerr=ystd, capsize=4, elinewidth=2, ecolor='k',
                marker='o', color='k', lw=0)
ax.plot(x, fit[0]*np.log10(x) + fit[1], color='red')
ax.set_xscale('log')
ax.set_xlabel('[Phenylalanine]/ M')
ax.set_ylabel('$k_{et}$/ $s^{-1}$')
ax.text(4e-5, 115, s)




fig, ax = plt.subplots(dpi=600)
ax.errorbar(concs, ket, yerr=ket_std, capsize=4, elinewidth=2, ecolor='k',
                marker='o', color='k', lw=0)
ax.set_xscale('log')
ax.set_xlabel('[Phenylalanine]/ M')
ax.set_ylabel('$k_{et}$/ $s^{-1}$')
ax.set_xticks([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2])