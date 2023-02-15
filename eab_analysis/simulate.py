import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
# plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')


# Initial params
# R1 = 100
# R2 = 10000
# Q1 = 1e-7
# Q2 = 1e-7
# n1 = 1
# n2 = 1


Rs = 100
R1 = 1000
R2 = 20000
Cd = 1e-7
C1 = 1e-7
C2 = 1e-7


f = np.logspace(0,np.log10(5000))

def CPE(f, params):
    '''
    Params:
        Q: CPE value
        n: CPE exponent (0 < n < 1)
    '''
    w = 2*np.pi*f
    Q = params['Q']
    n = params['n']
    return 1/(Q*(w*1j)**n)


# def Z(f, R1, R2, Q1, n1, Q2, n2):
#     '''
#     Params:
#         R1: Series resistance
#         R2: CT resistance
#         Q1: dl capacitance
#         n1: dl exponent (0 < n < 1)
#         Q2: Adsorption CPE
#         n2: Adsorption exponent (0 < n < 1)
#     '''
        
#     Ca = CPE(f, {'Q':Q2, 'n':n2})
#     Cdl = CPE(f, {'Q':Q1, 'n':n1})
    
#     Z = R1 + 1/(1/Cdl + 1/(R2+Ca))
    
#     return Z



def Z(f, Rs, R1, R2, Cdl, C1, C2):
    
    Cdl = CPE(f, {'Q':Cdl, 'n':1})
    C1 = CPE(f, {'Q':C1, 'n':1})
    C2 = CPE(f, {'Q':C2, 'n':1})
    
    z1 = R1 + C1
    z2 = R2 + C2
    
    Z = Rs + 1/((1/Cdl) + (1/z1) + (1/z2))
    
    return Z


def make_slider(pos, label, lim, valfmt=None):
    
    slider = Slider(
        ax      = plt.axes([0.65, pos, 0.25, 0.03]),
        label   = label,
        valmin  = lim[0],
        valmax  = lim[1],
        valinit = lim[2],
        valfmt = valfmt)
    
    return slider



    

# fig, ax = plt.subplots(figsize = (10,5), dpi=100)
# line,   = plt.plot(f, np.angle(Z(f, R1, R2, Q1, n1, Q2, n2), deg=True),
#                    color = colors[0])

# ax.set_xscale('log')
# ax.set_xlabel('Frequency/ Hz')
# ax.set_ylabel('Phase/ $\degree$')

# ax2 = ax.twinx()
# line2, = plt.plot(f, np.absolute(Z(f, R1, R2, Q1, n1, Q2, n2)), color=colors[1])
# ax2.set_ylabel('|Z|/ $\Omega$')
# ax2.set_yscale('log')

# plt.subplots_adjust(right=0.5)





def update(val):
    global ax, ax2
    
    # line.set_ydata(np.angle(Z(f, R1_slider.val,
    #                             R2_slider.val,
    #                             Q1_slider.val,
    #                             n1_slider.val,
    #                             Q2_slider.val,
    #                             n2_slider.val),
    #                deg=True))
    
    # line2.set_ydata(np.absolute(Z(f, R1_slider.val,
    #                             R2_slider.val,
    #                             Q1_slider.val,
    #                             n1_slider.val,
    #                             Q2_slider.val,
    #                             n2_slider.val)
    #                             ))
    
    line.set_ydata(np.angle(Z(f, *[slider.val for slider, _ in sliders]), deg=True))
    line2.set_ydata(np.absolute(Z(f, *[slider.val for slider, _ in sliders])))
    
    ax.relim()
    ax.autoscale_view()
    ax2.relim()
    ax2.autoscale_view()
    ax.figure.canvas.draw_idle()
    
    

# R1_slider = make_slider(0.7, 'Rs/ $\Omega$', [0, 5000, R1])
# R2_slider = make_slider(0.6, 'Rct/ $\Omega$', [100, 30000, R2])
# Q1_slider = make_slider(0.5, 'Cdl/ F', [1e-12, 1e-6, Q1], valfmt = '%.2E')
# Q2_slider = make_slider(0.4, 'Cad/ F', [1e-12, 1e-6, Q2], valfmt = '%.2E')
# n1_slider = make_slider(0.3, 'n_dl', [0, 1, n1])
# n2_slider = make_slider(0.2, 'n_ad', [0, 1, n2])

Rs_slider = None
R1_slider = None
R2_slider = None
Cd_slider = None
C1_slider = None
C2_slider = None
Q1_slider = None
Q2_slider = None
n1_slider = None
n2_slider = None

# sliders = [(R1_slider, (0.7, 'Rs/ $\Omega$', [0, 5000, R1])),
#            (R2_slider, (0.6, 'Rct/ $\Omega$', [100, 30000, R2])),
#            (Q1_slider, (0.5, 'Cdl/ F', [1e-12, 1e-6, Q1], '%.2E')),
#            (n1_slider, (0.4, 'Cad/ F', [1e-12, 1e-6, Q2], '%.2E')),
#            (Q2_slider, (0.3, 'n_dl', [0, 1, n1])),
#            (n2_slider, (0.2, 'n_ad', [0, 1, n2]))]

ss = [(Rs_slider, (0.7, 'Rs/ $\Omega$', [0, 5000, Rs])),
           (R1_slider, (0.6, 'R1/ $\Omega$', [0, 50000, R1])),
           (R2_slider, (0.5, 'R2/ $\Omega$', [0, 50000, R2])),
           (Cd_slider, (0.4, 'Cdl/ F', [1e-12, 1e-6, Cd], '%.2E')),
           (C1_slider, (0.3, 'C1/ F', [1e-12, 1e-6, C1], '%.2E')),
           (C2_slider, (0.2, 'C2/ F', [1e-12, 1e-6, C2], '%.2E'))]




fig, ax = plt.subplots(figsize = (10,5), dpi=100)

sliders = []

for slider, args in ss:
    slider = make_slider(*args)
    slider.on_changed(update)
    
    sliders.append((slider, args))



plt.sca(ax)

line,   = plt.plot(f, np.angle(Z(f, *[slider.val for slider, _ in sliders]), 
                               deg=True),
                   color = colors[0])

ax.set_xscale('log')
ax.set_xlabel('Frequency/ Hz')
ax.set_ylabel('Phase/ $\degree$')

ax2 = ax.twinx()
line2, = plt.plot(f, np.absolute(Z(f, *[slider.val for slider, _ in sliders])), 
                  color=colors[1])
ax2.set_ylabel('|Z|/ $\Omega$')
ax2.set_yscale('log')

plt.subplots_adjust(right=0.5)


plt.show()




