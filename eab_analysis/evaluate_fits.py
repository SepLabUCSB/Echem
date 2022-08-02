import os
import numpy as np
import matplotlib.pyplot as plt


data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-04-13\6 hr with dosing -330mV 1154am'

folders = []
for file in os.listdir(data_dir):
    if len(file.split('.')) > 1:
        continue
    else:
        folders.append(file)
        
folders.sort(key=lambda s:int(s))

re_deltas = []
im_deltas = []
for folder in folders:
    folder = os.path.join(data_dir, folder)
    for file in os.listdir(folder):
        if file.endswith('_fit.txt'):
            fit_file = file
            base_file = f'{file[:4]}s.txt'
            
            f, re, im         = np.loadtxt(os.path.join(folder, base_file), 
                                           unpack=True, skiprows=1)
            f, re_fit, im_fit = np.loadtxt(os.path.join(folder, fit_file), 
                                           unpack=True)
            
            re     = re[:21]
            im     = im[:21]
            re_fit = re_fit[:21]
            im_fit = im_fit[:21]
            
            re_deltas.append((re_fit - re)/re)
            im_deltas.append((im_fit - im)/im)
 

rd = np.average(np.array(re_deltas), axis=0)
imd = np.average(np.array(im_deltas), axis=0)
            
fig, ax = plt.subplots()
ax.plot(f[:21], rd, label='Re')
ax.plot(f[:21], imd, label='Im')
ax.set_xscale('log')
ax.set_xlabel('Frequency/ Hz')
ax.set_ylabel('Relative fit error')
ax.show_legend()