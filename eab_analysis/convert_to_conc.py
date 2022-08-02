import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
file = 'C:/Users/BRoehrich/Desktop/vanco data/meisp_fits.xlsx'

# def f(k, n, Keq, A, B):
#     a = (k-A)/B
#     return ((a*Keq)/(1-a))**(1.0/n)


def f(k, n, KD, A, B):
    return ((KD*(k-A))/(A+B-k))**(1/n)


def hill_equation(C, n, Keq, A, B):
    # n = 1
    return A + B * ((C**n)/(Keq + C**n))


n  = 1.14008096e+00
KD = 3.67633833e-05
A  = 2.15598378e+01
B  = 3.12810297e+01
# n, KD, A, B = [1.00000000e+00, 1.41934290e-04, 2.13062958e+01, 3.32865076e+01]


# t,Rs,Rct,Cdl,Cad,ket,_,_,_,_ = np.loadtxt(file, unpack=True, skiprows=1)
df = pd.read_excel(file, skiprows=1, names=('t','Rs','Rct','Cad','phi','Cdl', 'ket',
                                            'a', 'b', 'c', 'd', 'e'))

t = df['t'].to_numpy()
Rs = df['Rs'].to_numpy()
Rct = df['Rct'].to_numpy()
Cdl = df['Cdl'].to_numpy()
Cad = df['Cad'].to_numpy()
ket = df['ket'].to_numpy()
conc = f(ket, n, KD, A, B)


fig, ax = plt.subplots(figsize=(5,5), dpi=300)
ax.plot(t/60, conc/1e-6, '.')
ax.set_xlabel('Time/ min')
ax.set_ylabel('[Vancomycin]/ $\mu$M')


# fig, ax = plt.subplots(figsize=(5,5), dpi=100)
# ax.plot(t/60, ket, '.')
# ax.set_xlabel('Time/ min')
# ax.set_ylabel('$k_{et}$/ $s^{-1}$')