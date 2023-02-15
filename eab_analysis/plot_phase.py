import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


# data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-09-01\-340mV vanco'
data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-07-27 vanco invivo'
# data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-04-05'

def extract_phase_data(data_dir, probe=''):
    
    data = []
    nums = []
    
    # time_file = data_dir + '/0000_time_list.txt'
    # times = np.loadtxt(time_file)   
    
    for file in os.listdir(data_dir) :
        if file.startswith(probe) and file.endswith('.txt'):
            if ('currents' in file or 'Metadata' in file or 'time_list' in file):
                continue
            # num = int(file.strip('.txt').split('_')[-1])
            num = int(file.strip('s.txt'))
            f, re, im = np.loadtxt(os.path.join(data_dir, file), unpack=True,
                                   skiprows=1)
            
            Z = re + 1j*im
            
            data.append((num, Z))
            nums.append(num)
            
    
    data.sort(key=lambda x:x[0])
    nums.sort()
    
    _, data = zip(*data)
    
    # times = times[tuple([nums])] # np weirdness, get times from indices
    
    times = nums
    
    return times, f, data



def plot_phase(freq, times, freqs, data):
    
    ind = np.where(freqs == freq)[0][0]
    
    p = [spec[ind] for spec in data]
    p = np.angle(p, deg=True)
    
    fig, ax = plt.subplots()
    ax.plot(times, p)
    ax.set_xlabel('Time/ s')
    ax.set_ylabel(f'Phase @ {freq:.0f} Hz')
    
    return p


def plot_phase_ratio(times, data, freqs, f1, f2):
    
    ind1 = np.where(freqs == f1)[0][0]
    ind2 = np.where(freqs == f2)[0][0]
    
    p1 = [spec[ind1] for spec in data]
    p1 = np.angle(p1, deg=True)
    
    p2 = [spec[ind2] for spec in data]
    p2 = np.angle(p2, deg=True)
    
    fig, ax = plt.subplots()
    # ax.plot(times, p1/p1[0])
    # ax.plot(times, p2/p2[0])
    ax.plot(times, (p1/p1[0])/(p2/p2[0]))
    ax.set_xlabel('Time/ s')
    ax.set_ylabel(f'Phase({f2:.0f})/ Phase({f1:.0f})')
    
    return 


def plot_Bodes(times, freqs, data):
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(data)))
    
    i = 0
    fig, ax = plt.subplots()
    for d in data:
        
        if i%100 == 0:
            ax.plot(freqs, np.angle(d, deg=True), color=colors[int(i)])
        
        
        i += 1
    ax.set_xscale('log')
    ax.set_xlabel('Frequency/ Hz')
    ax.set_ylabel('Phase/ $\degree$')


times, freqs, data = extract_phase_data(data_dir)
plot_Bodes(times, freqs, data)
# for freq in freqs:
#     plot_phase(freq, times, freqs, data)
    
    
    
    
    
    
    
    
    
    
    
    