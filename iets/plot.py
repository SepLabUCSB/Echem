import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os


# df = pd.read_csv('C:/Users/BRoehrich/Desktop/iets/6MHA-Au 1mM FcMeOH 100mM KCl PtCE SCE17/CVs/CV 6MHA Au 1mM FcMeOH pre IETS.txt',
#                  names=('a', 'b', 'i', 'v', '1', '2', '3', '4', '5', '6'),
#                  sep=';', skiprows=1)

# i = df['i'].to_numpy()
# v = df['v'].to_numpy()


folder = r'C:\Users\BRoehrich\Desktop\iets\6-MHOH-Au Ruhex 100 scans'
# monofolder = r'C:\Users\BRoehrich\Desktop\iets\6MHA-Au 1mM FcMeOH 100mM KCl PtCE SCE17'

i = 0
pol = []
for file in os.listdir(folder):
    
    try:
        i = int(file.strip('.txt'))
    except:
        continue
    
    if i < len(os.listdir(folder))-1:
        v, I = np.loadtxt(os.path.join(folder,file), skiprows=1, unpack=True,
                          delimiter=',')
        pol.append((v,I))
    
    # i += 1

# i = 0
# mono = []

# for file in os.listdir(monofolder):
    
#     try:
#         i = int(file.strip('.txt'))
#     except:
#         continue
    
#     if i < len(os.listdir(monofolder))-1:
#         v, I = np.loadtxt(os.path.join(monofolder,file), skiprows=1, unpack=True,
#                           delimiter=',')
#         mono.append((v,I))
    
#     # i += 1


#%%
    
v, i = zip(*pol)
# v = list(v)
# i = list(i)
# v.pop(-3)
# i.pop(-3)
v = np.average(v, axis=0)
i = np.average(i, axis=0)


# v2, i2 = zip(*mono)
# v2 = np.average(v2, axis=0)
# i2 = np.average(i2, axis=0)

# diff = i2 - i

# v = v[:int(len(v)/2)]
# i = i[:int(len(i)/2)]

# v = np.mean(np.reshape(v, (-1,5)), axis=1)
# i = np.mean(np.reshape(i, (-1,5)), axis=1)



fig, ax = plt.subplots(figsize=(5,5), dpi=100)
ax.plot(v, i, '.', label='Polished')
# ax.plot(v2, i2, '.', label='Monolayer')
# ax.plot(v, diff)
ax.set_xlabel('E/ V vs SCE')
ax.set_ylabel('d2I/dV2')
# ax.legend()
plt.show()
# 

# fig, ax = plt.subplots(figsize=(5,5), dpi=100)
# ax.plot(v[1:], np.diff(i, n=1), '.')
# ax.set_xlabel('E/ V vs SCE')
# ax.set_ylabel('dI')
# # plt.show()

# fig3, ax = plt.subplots(figsize=(5,5), dpi=100)
# ax.plot(v[1:], np.diff(i, n=1), '.-')
# ax.set_xlabel('E/ V vs SCE')
# ax.set_ylabel('d2I/dV2')
# plt.show()

#%%
# df = pd.read_csv('C:/Users/BRoehrich/Desktop/iets/Au ORR 100mM KCl vs SCE 100scans wider window/Au ORR CV.txt',
#                  names=('a', 'b', 'i', 'v', '1', '2', '3', '4', '5', '6'),
#                  sep=';', skiprows=1)

# i = df['i'].to_numpy()
# v = df['v'].to_numpy()


# a = 186
# b = -1228
# v = v[a:b]
# i = i[a:b]


# fig, ax = plt.subplots()
# ax.plot(v, 1e6*i, label='I')
# ax.plot(v[1:], 13e6*np.diff(i, n=1) - 10, label='dI/dV')
# ax.plot(v[2:], 5e6*np.diff(i, n=2) - 15, label='$d^{2}I/dV^{2}$')
# ax.set_xlabel('E/ V vs SCE')
# ax.set_ylabel('I or dI/dV or $d^{2}I/dV^{2}$')
# ax.set_yticks([])
# ax.legend()