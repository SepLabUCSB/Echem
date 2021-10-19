# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 09:44:33 2021

@author: Condor
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import count
import warnings
warnings.filterwarnings("ignore")
#plt.rcParams['figure.dpi'] = 300

class Spike:
    _ids = count(0)
    
    def __init__(self, peak, time):
        self.id = next(self._ids) #count spike number
        self.peak = peak
        self.time = time
        
    
    def peak_integrate(self, data, dx):
        points_duration = self.right - self.left
        time_duration = points_duration*dx
        # print(points_duration, 'points,', time_duration, 's')
        baseline_area = (data[self.left]+data[self.right])*time_duration*0.5
        integral = np.trapz(data.loc[self.left:self.right], dx=dx)
        self.integral = integral - baseline_area


'''
https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/43512887#43512887
'''

def thresholding_algo(y, lag, threshold, influence):
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] < avgFilter[i-1]:                                   # spikes be negative boi
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))        


#df = pd.read_fwf('C:/Users/BRoehrich/Desktop/lsv4.txt', names=('time', 'current'), widths=[23,23])
df = pd.read_csv('D:/_DataProcessing/In use/1..txt', delimiter = '\t', skiprows = 1, names = ('time','potential','current','range'))
#df = df[1:]

df['current'] = df['current'].astype(float)
df['time'] = df['time'].astype(float)
df = df[df['time']>1]
df = df[df['time'] > 50]
df['current'] = df['current']*1e-3
# df['current'] = df['current'] + 1
current = df['current'].to_numpy()
time = df['time'].to_numpy()
dx = time[1] - time[0]

#%%

out = thresholding_algo(current, 500, 3, 0.01)
points = out['signals']
avg = out['avgFilter']
df['spike'] = points
df['avg'] = avg

spikes = {}

## Delete extraneous points (tags multiple points in the same peak, reduce to 1 per peak)
for i in range(df.index[0], df.index[0]+len(df)):
    if df.loc[(i, 'spike')] == 1:
        if i-1 not in spikes and i-2 not in spikes and i-3 not in spikes and i-4 not in spikes and i-5 not in spikes and i-6 not in spikes and i-7 not in spikes and i-8 not in spikes and i-9 not in spikes and i-10 not in spikes and i-11 not in spikes and i-12 not in spikes and i-13 not in spikes and i-14 not in spikes and i-15 not in spikes:
            time_i = df.loc[(i, 'time')]
            spikes[i] = Spike(i, time_i)
    
hist = []   
for i in spikes:
    ## Set left bound for integration
    ## Boundary where current falls to less than the running average from thresholding_algo
    n_left = 0
    c = abs(df.loc[(i+n_left, 'current')])
    while abs(c) > abs(df.loc[i, 'avg']):
        n_left = n_left-1
        c = df.loc[(i+n_left, 'current')]
        
    
    n_right = 0
    c = abs(df.loc[(i+n_right, 'current')])
    while c > abs(df.loc[i, 'avg']):
        n_right = n_right + 1
        try:
            c = df.loc[(i+n_right, 'current')]
        except:
            print('Ran into end of dataset!')
            n_right = n_right-1
            break
        # print(c, abs(df.loc[i, 'avg']))
    
    spikes[i].left = int(i + n_left)
    spikes[i].right = int(i + n_right)
    spikes[i].peak_integrate(df['current'], dx=dx)
    hist.append(1e15*spikes[i].integral)

print('Integrated %s spikes.' %len(spikes))

#%%
plt.figure()
plt.plot(time, current, '.-')
plt.xlabel('Time/ s')
plt.ylabel('Current/ A')
for i in spikes:
    plt.plot(spikes[i].time, df.loc[(i, 'current')], 'xr')
    #plt.annotate('A spike!', xy = (spikes[i].time, df.loc[(i, 'current')]), arrowprops = {'color':'orange'})
    #print(i)
plt.plot(time, out['avgFilter'], '-r', label = 'avgFilter')
plt.plot(time, out['avgFilter'] + out['stdFilter'], '-g', label = '+std')
plt.plot(time, out['avgFilter'] - out['stdFilter'], '-g', label = '-std')
plt.xlim(90,100)
plt.ylim(-4.5e-10, -4e-10)
plt.legend()
plt.show()




plt.hist(hist, bins=np.arange(-10000,-200,500), rwidth=0.8)
plt.xlabel('Charge/ fC')
plt.ylabel('Count')