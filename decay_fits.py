import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from scipy import optimize, signal
import os
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
data_folder = r'C:\Users\BRoehrich\Desktop\Miguel data'


correct_baselines = False
start_after = 100 # cut off first (n) seconds
delay = 10 # Points after "fast" spike to skip fitting on

apply_filter     = True
filter_freq      = 25

fit = 'monoexponential'
# fit = 'biexponential'
# fit = 'linear'



# Define fit function and parameters    
if fit == 'biexponential':
    params=['a/ A',
            'b/ s-1',
            'c/ A',
            'd/ s-1',
            'e/ A']
    
    def exp_func(x, a,b,c,d,e):
        return a * np.exp(-b * x) + c * np.exp(-d * x) + e


elif fit == 'monoexponential':
    params = ['a/ A',
            'b/ s-1',
            'c/ A']
    
    def exp_func(x, a,b,c):
        return a * np.exp(-b * x) + c


elif fit == 'linear':
    params = ['a/ A s-1', 
              'b/ A']
    
    def exp_func(x, a, b):
        return a*x + b



def main(folder, params):
    # Analyze all files in folder and save together   
    d = {param: [] for param in params}
    d['integral/ C'] = []
    d['Reduced Chi^2'] = []
    d['time/ s'] = []
    d['File'] = []
    
    
    for file in os.listdir(folder):    
        if file.endswith('.txt'):
            # fits = list of lists i.e. [[*param1], [*param2], ...]
            # integrals: spike areas. Unit: C
            fits, integrals, chi_sqs, points = analyze_file(folder + '/' + file,
                                                    correct_baselines, 
                                                    apply_filter, delay,
                                                    start_after)
            
            for i in range(len(params)):
                for x in fits[i]:
                    d[params[i]].append(x)
            for x in integrals:
                d['integral/ C'].append(x)
            for x in chi_sqs:
                d['Reduced Chi^2'].append(x)
            for x in points:
                d['time/ s'].append(x)
            
            d['File'].append(file)
            d['File'].append(f'Fit: {fit}')
            d['File'].append(f'Baseline correct: {correct_baselines}')
            d['File'].append(f'Delay: {delay} pts')
            
            for key in d:
                while len(d[key]) < max([len(d[key]) for key in d]):
                    d[key].append(' ')


    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in d.items() ]))

    
    writer = pd.ExcelWriter(folder + '/output.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {folder}/output.xlsx')
    
    return d, df






# Function for reading data
def extract_data(file, start_after):
    df = pd.read_csv(file, names=('t', 'v', 'i'), skiprows=1, sep='\t')
    df = df[df['t'] > start_after]
    
    t = np.array(df['t'])
    i = np.array(df['i'])/1000
    
    # Calculate sampling rate
    sara = np.average([t[i] - t[i-1] for i in range(1, len(t))])
    
    return t, i, sara




# Filtering functions
def lowpass(y, t, cutoff):
    
    fs = 1/np.mean([t[i] - t[i-1] for i in range(1, len(t))])
    
    fc = cutoff/(fs/2)
    
    try:
        b, a = signal.butter(8, fc)
        filt_y = signal.filtfilt(b, a, y, padlen=150)
        return filt_y
    
    except ValueError:
        print('Bad filter_freq, not filtering.')
        print(f'Valid filter_freq from 0 < f < {fs/2}')
        return y
    






def thresholding_algo(y, lag, threshold, influence):
    '''
    Taken from 
    https://stackoverflow.com/questions/22583391/peak-signal-detection-in-
    realtime-timeseries-data/43512887#43512887
    
    
    Returns:
        signals: array of [-1, 0, 1]. 0 = no spike at this point,
                 +- 1 = positive/ negative spike at this point
        avgFilter: array, filtered data
        stdFilter: array, standard devation of filtered data
    '''
    
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
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
    
    
    return signals, np.array(avgFilter), np.array(stdFilter)



def refine_peaks(y, signals):
    '''
    Finds local maxima to refine peak locations

    Parameters
    ----------
    y : array
        Raw data.
    signals : array
        Output from thresholding_algo, array of [-1, 0, 1]. 
        Must be same length as y.

    Returns
    -------
    idxs : list
        List of indices corresponding to spike locations.

    '''
    
    # Determine index of this peak and next peak
    idxs = list(np.where(signals != 0)[0])
    
    # Remove double labelled spikes
    for idx in idxs:
        for i in range(idx-10, idx+10):
            if i in idxs and i != idx:
                print(f'Removing {i} because of duplicate')
                idxs.remove(i)
        
        
    # Refine peak location            
    for idx in idxs[:]:
        if any(abs(y[idx-10:idx+10]) > abs(y[idx])):
            
            i = np.where(abs(y) ==
                          max(abs(y[idx-10:idx+10])))[0][0]
            
            print(f'Moving {idx} to {i}')
            idxs.remove(idx)
            idxs.append(i)
            idxs.sort()
            
    return idxs



def integrate_spikes(t, y, idxs, avg):
    '''
    Parameters
    ----------
    t : array
        Time.
    y : array
        Raw data.
    idxs : list
        List of spike indices.
    avg : array
        Filtered current data from thresholding_algo().

    Returns
    -------
    integrals : list
        List of peak integrals. Unit Amperes.

    '''
    
    integrals = []
    
    
    for idx in idxs:
        
        # Determine left and right bounds for integration.
        # Defaults to peak location +- 1 point
        # Try to refine bounds by looking for where the data crosses
        # the filtered "avg" line (prev used to detect where spikes are),
        # if there are crossings within +- 3 points, use those as the 
        # bounds instead.
        
        try:
            left_bound = idx - 3 + np.where(abs(y[idx-3:idx]) <
                                      abs(avg[idx-3:idx]))[0][-1]              
        except:
            left_bound = idx-1
        
        try:
            right_bound = idx + np.where(abs(y[idx:idx+3]) <
                                      abs(avg[idx:idx+3]))[0][0]              
        except:
            right_bound = idx+1
        
        
        bounds = [left_bound, right_bound]
        
        y1 = y[bounds[0]]
        y2 = y[bounds[1]]
        dt = t[bounds[1]] - t[bounds[0]]
        
        baseline_area = (1/2) * (y1 + y2) * dt
        
        total_area = np.trapz(y[bounds[0]:bounds[1]],
                               x = t[bounds[0]:bounds[1]])
                
        
        integral = total_area - baseline_area
        integrals.append(integral)
        
                
        
        # fig = plt.figure(figsize=(5,5), dpi=100)
        # plt.plot(t[idx-500:idx+500], y[idx-500:idx+500])
        # plt.plot(t[idx], y[idx], 'ro')
        # plt.plot(t[idx-500:idx+500], avg[idx-500:idx+500])
        # plt.plot(t[bounds], y[bounds], 'o', color = 'orange')
                
        # print(integral)
    
    return integrals



def fit_peaks(t, y, idxs, sara, ax=None, t_max=10, delay = 10,
              baseline_correct=True, app_filter=False): 
    '''
    Fit exponential decay betweek this peak and next peak,
    or the next n seconds

    Parameters
    ----------
    t : array
        Time.
    y : array
        Raw data.
    idxs : list
        List of spike indices.
    sara : float
        Sampling rate from extract_data().
    ax : plt.Axes, optional
        ax to plot exponential fits on top of.
    t_max : int, optional
        Maximum time to use for fitting each spike. The default is 10.
    delay : int, optional
        Number of points after fast spike to skip for curve fitting

    Returns
    -------
    fits : list
        List of fitted [a, b, c] parameters for each spike.
    '''
    
    fits = []
    chi_sqs = []
    pops = []
    
    if app_filter:
        y = lowpass(y, t, filter_freq)
    
    for i, _ in enumerate(idxs):
        
        # Get this index, either use next index or + (delay) pts        
        this_idx = idxs[i] + delay
        try:
            next_idx = min(idxs[i+1] - 2, this_idx + int(round(t_max/sara)))
        except: # Last point
            next_idx = min(len(y), this_idx + int(round(t_max/sara)))
        
            
        # Subtract out baseline
        # Either use previous 500 pts, or last index
        
        if i > 0:
            last_idx = max(idxs[i-1] + delay, this_idx - 500)
        elif i == 0: # first spike
            last_idx = max(0, this_idx - 500)
        
        if (this_idx-delay - last_idx) < 100:
            print(f'Not enough points to fit baseline {t[this_idx-delay]} s.')
            pops.append(this_idx - delay)
            continue
        
        m, b = np.polyfit(t[last_idx:this_idx-delay-5], 
                  y[last_idx:this_idx-delay-5], deg=1)
        
        # Subtract the baseline
        if baseline_correct:
            baseline = m*t + b
        else:
            baseline = 0*t
        
        if (next_idx - this_idx) < 50:
            print(f'Not enough points to fit {t[this_idx - delay]} s.')
            pops.append(this_idx - delay)
            continue
            
        ts = t[this_idx:next_idx] - t[this_idx]
        data = y[this_idx:next_idx] - baseline[this_idx:next_idx]
        
        if fit == 'biexponential':
            bounds=(-np.inf, [1., np.inf, np.inf, np.inf, np.inf])
        elif fit == 'monoexponential':
            bounds=(-np.inf, [1., np.inf, np.inf])
            
        popt, pcov = optimize.curve_fit(exp_func, ts, data, 
                            maxfev=100000,
                            # bounds=bounds
                            )
        
        print(popt)        
        fits.append(popt)
        
        
        # Calculate chi^2
        fit_y = exp_func(ts, *popt) # Fits
        residuals = abs((data - fit_y)/fit_y)
        chi_sq = np.sum(residuals**2)/len(residuals)
        chi_sqs.append(chi_sq)
        
        # fig, ax2 = plt.subplots(figsize=(5,5), dpi=100)
        # ax2.plot(ts, data, '.')
        # ax2.plot(ts, fit_y)
        
        
        ax.plot(ts+t[this_idx], fit_y + baseline[this_idx:next_idx], 'gold')
        if baseline_correct:
            ax.plot(t[last_idx:next_idx], baseline[last_idx:next_idx], 'r--')
        
    
           
        
    return fits, chi_sqs, pops



def analyze_file(file, baseline_correct, app_filter, delay, start_after):
    '''
    Main function to call to detect and fit spikes
    '''
    
    # Get data
    t, i, sara = extract_data(file, start_after)
    
    # Find peaks
    signals, avgFilter, stdFilter = thresholding_algo(i, lag=50, threshold=5,
                                                        influence = 0.6)
    # Refine peak locations
    idxs = refine_peaks(i, signals)
    
    
    # Make figure
    if 'inline' in matplotlib.get_backend():
        fig, ax = plt.subplots(figsize=(5,5), dpi=300)
    else:
        fig, ax = plt.subplots(figsize=(5,5), dpi=100)
    
    if app_filter:
        ax.plot(t, lowpass(i,t,filter_freq))
    else:
        ax.plot(t, i)
    ax.plot(t[idxs], 1.01*i[idxs], 'ro')
    ax.set_xlabel('Time/ s')
    ax.set_ylabel('Current/ A')
    plt.tight_layout()
    
    
    # Fit peaks. Optionally plot fits onto same ax    
    fits, chi_sqs, pops = fit_peaks(t, i, idxs, sara, ax, 
                              baseline_correct=baseline_correct,
                              app_filter = app_filter, delay=delay)
    
    for pt in pops:
        # Remove spike indices that we couldn't fit
        idxs.pop(np.where(idxs==pt)[0][0])
    
    # Calculate spike integrals
    integrals = integrate_spikes(t, i, idxs, avgFilter)
    
    # Get time points
    points = t[idxs]
    
    # Transpose data for saving
    fits = np.array(fits).T
    a, b, c = fits[0], fits[1], fits[2]
    
    
    return fits, integrals, chi_sqs, points







d, df = main(data_folder, params)         
        
    
