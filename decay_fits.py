import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from scipy import optimize, signal
import os
import warnings
plt.ion()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# plt.style.use('C:/Users/orozc/Google Drive (miguelorozco@ucsb.edu)/Research/Spyder/scientific.mplstyle')
# data_folder = r'C:/Users/orozc/Google Drive (miguelorozco@ucsb.edu)/Research/Spyder/Run'
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
data_folder = r'C:\Users\BRoehrich\Desktop\Miguel data'



correct_baselines = False
start_after = 100 # cut off first (n) seconds
min_s_to_fit = 10  # Requires n seconds of data to accept the fit
delay = 10 # Points after "fast" spike to skip fitting on
thresh = 0.01 # Used to determine acceptable baseline "flatness"
               # Smaller = more picky, need flatter baseline to accept spike

apply_filter     = False
filter_freq      = 25

# fit = 'monoexponential'
fit = 'monoexp-linear'
# fit = 'biexponential'
# fit = 'linear'



''' 
TO-DO
    
Interactive picking/deleting points

'''


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
    
elif fit == 'monoexp-linear':
    params = ['a/ A',
            'b/ s-1',
            'm/ A s-1',
            'c/ A']
    
    def exp_func(x, a, b, m, c):
        return a*np.exp(-b * x) + m*x + c



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




class InteractivePicker:
    
    global colors
    
    def __init__(self, file, fig, ax, params, exp_func, plot=True):
        self.file = file
        self.fig  = fig
        self.ax   = ax
        self.plot = plot
        self.params = params
        self.exp_func = exp_func
        
        self.fits   = []
        self.points = []
        
        self.line        = None
        self.fitline     = None
        self.drawnpoints = None
        
        self.save_fits  = []
        self.integrals  = []
        self.chi_sqs    = []
        self.points     = []
        self.ms         = []
        
    
    def __call__(self, event):
        if not self.plot:
            return
        # Set xtol based on current axis limits 
        upper = self.ax.get_xlim()[1]
        lower = self.ax.get_xlim()[0]
        diff = upper-lower
        xtol = 0.01*diff
        
        
        if event.inaxes and self.fig.canvas.manager.toolbar.mode == "":
            # print(xtol)
            clickX = event.xdata
            if (self.ax is None) or (self.ax is event.inaxes):
                # Prioritizing removing marked point
                for x in self.t[self.idxs]:
                    if (clickX-xtol < x < clickX+xtol):
                        # Remove point
                        print(f'Removing closest point: {x}')
                        idx = np.where(self.t == x)[0][0]
                        i   = np.where(self.idxs == idx)[0][0]
                        self.idxs.remove(idx)
                        
                        t = [(i, x, y) for (i, x, y) in self.fits if i == idx][0]
                        self.fits.remove(t)
                        
                        self.draw_fits(xlim=self.ax.get_xlim(),
                                       ylim=self.ax.get_ylim())
                        
                        
                        # Remove from output lists
                        self.save_fits  = np.delete(self.save_fits, i, axis=1)
                        self.areas      = np.delete(self.areas, i)
                        self.integrals  = np.delete(self.integrals, i)
                        self.chi_sqs    = np.delete(self.chi_sqs, i)
                        self.points     = np.delete(self.points, i)
                        self.ms         = np.delete(self.ms, i)
                        
                    else:
                        continue
    
    
    def keypress_event_handler(self, event):
        # Pressing arrow keys pans around graph
        key = event.key
        if key in ['left', 'right', 'up', 'down']:
            shift_axes(self.ax, key)
        self.ax.figure.canvas.draw_idle()
    
    
    def scroll_event_handler(self, event):
        # Scrolling mousewheel zooms in/out
        key = event.button
        x, y = event.xdata, event.ydata
        if key in ['down', 'up']:
            zoom_axes(self.ax, key, (x,y))
        self.ax.figure.canvas.draw_idle() 
                
                    
    def _reset(self, event):
        file     = self.file
        fig      = self.fig
        ax       = self.ax
        params   = self.params
        exp_func = self.exp_func
        
        self.line.remove()
        self.fitline.remove()
        self.drawnpoints.remove()
        self.ax.clear()
        
        self.fig.canvas.draw()
        
        for attr, _ in list(self.__dict__.items()):
            delattr(self, attr)
        
        self.__init__(file, fig, ax, params, exp_func)
        self.fig.canvas.mpl_connect('button_press_event', self)
        
        fits, integrals, chi_sqs, points, ms = self.analyze_file(
                                                    correct_baselines, apply_filter, 
                                                    delay, start_after, self.params,
                                                    self.exp_func, thresh
                                                    )
        
        
        return
    
    
    
    def draw_fits(self, xlim=[], ylim=[]):
        if not self.plot:
            return
        # Plot fit curves
        if self.fitline:
            self.fitline.remove()
        
        self.fitline, = self.ax.plot([0], [0], 'o', color='gold',
                                     ms=1)
        
        xs = np.array([])
        ys = np.array([])
        
        for (idx, x_data, y_data) in self.fits:
            xs = np.hstack([xs, x_data])
            ys = np.hstack([ys, y_data])
        
        self.fitline.set_data(xs, ys)
        
        # Plot point labels
        if self.drawnpoints:
            self.drawnpoints.remove()
        self.drawnpoints, = self.ax.plot([0], [0], 'ro')
        
        xs = []
        ys = []
        for idx in self.idxs:
            xs.append(self.t[idx])
            ys.append(0.95*self.i[idx])
        
        self.drawnpoints.set_data(xs, ys)        
        
        self.fitline.figure.canvas.draw()
        
        if len(xlim) != 0:
            self.ax.set_xlim(xlim)
        if len(ylim) != 0:
            self.ax.set_ylim(ylim)
            
        return
    
    
    # Function for reading data
    def extract_data(self, file, start_after):
        df = pd.read_csv(file, names=('t', 'v', 'i'), skiprows=1, sep='\t')
        df = df[df['t'] > start_after]
        
        t = np.array(df['t'])
        i = np.array(df['i'])/1000
        
        # Calculate sampling rate
        sara = np.average([t[i] - t[i-1] for i in range(1, len(t))])
        
        self.t = t
        self.i = i
        
        return t, i, sara
          
    
    def refine_peaks(self, y, signals):
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
                    #print(f'Removing {i} because of duplicate')
                    idxs.remove(i)
            
            
        # Refine peak location            
        for idx in idxs[:]:
            if any(abs(y[idx-10:idx+10]) > abs(y[idx])):
                
                i = np.where(abs(y) ==
                              max(abs(y[idx-10:idx+10])))[0][0]
                
                #print(f'Moving {idx} to {i}')
                idxs.remove(idx)
                idxs.append(i)
                idxs.sort()
                
        return idxs
    
    
    
    def integrate_spikes(self, t, y, idxs, avg):
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
        pops      = []
        
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
            
            if integral < 0:
                pops.append(idx)
                print(f'Negative integral: {idx}')
            
            else:
                integrals.append(integral)
            
                    
            
            # fig = plt.figure(figsize=(5,5), dpi=100)
            # plt.plot(t[idx-500:idx+500], y[idx-500:idx+500])
            # plt.plot(t[idx], y[idx], 'ro')
            # plt.plot(t[idx-500:idx+500], avg[idx-500:idx+500])
            # plt.plot(t[bounds], y[bounds], 'o', color = 'orange')
                    
            # print(integral)
        
        return integrals, pops
    
    
    
    def fit_peaks(self, t, y, idxs, sara, ax=None, t_max=60, t_min=min_s_to_fit, 
                  delay = 10, baseline_correct=True, app_filter=False, 
                  thresh= 0.005): 
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
        t_min: int, optional
            Minimum time necessary to choose a spike for fitting. Default 1 s
        delay : int, optional
            Number of points after fast spike to skip for curve fitting
        thresh: float, optional
            Used to determine acceptable baseline "flatness"
    
        Returns
        -------
        fits : list
            List of fitted [a, b, c] parameters for each spike.
        '''
        
        fits    = []
        chi_sqs = []
        integrals = []
        pops    = []
        ms      = []
        
        if app_filter:
            y = lowpass(y, t, filter_freq)
        
        for i, _ in enumerate(idxs):
            
            # Get this index, either use next index or + (delay) pts        
            this_idx = idxs[i] + delay
            try:
                next_idx = min(idxs[i+1] - 2, this_idx + int(round(t_max/sara)))
            except: # Last point
                next_idx = min(len(y)-1, this_idx + int(round(t_max/sara)))
            
                
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
            
            
            
            # Calculate the baseline
            m, b = np.polyfit(t[last_idx:this_idx-delay-5], 
                  y[last_idx:this_idx-delay-5], deg=1)
            
            
            if abs(m/y[last_idx]) > thresh:
                print(f'Baseline not flat enough {t[this_idx-delay]} s.')
                pops.append(this_idx - delay)
                continue
            
            if baseline_correct:    
                baseline = m*t + b
            else:
                baseline = 0*t + 0
                
                
            
            
            if (next_idx - this_idx) < int(round(t_min/sara)):
                print(f'Not enough points to fit {t[this_idx - delay]} s.')
                pops.append(this_idx - delay)
                continue
            
            
            ### Do fitting ###
            ts = t[this_idx:next_idx] - t[this_idx]
            data = y[this_idx:next_idx] - baseline[this_idx:next_idx]
            
            if fit == 'biexponential':
                bounds=(-np.inf, [1., np.inf, np.inf, np.inf, np.inf])
            elif fit == 'monoexponential':
                bounds=(-np.inf, [1., np.inf, np.inf])
            
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    popt, pcov = optimize.curve_fit(exp_func, ts, data, 
                                    maxfev=100000,
                                    # bounds=bounds
                                    )
            except: 
                popt = [0,0,0,0,0]
            
                       
            
            # Calculate chi^2
            fit_y = exp_func(ts, *popt) # Fits
            residuals = abs((data - fit_y)/fit_y)
            chi_sq = np.sum(residuals**2)/len(residuals)
            
            
            ### Calculate integral (excludes initial sharp spike)      
            # area to draw
            verts = [(t[this_idx], fit_y[-1]),
                     *zip(t[this_idx:next_idx], data),
                     (t[next_idx], fit_y[-1])]
            poly = Polygon(verts, color= colors[0], alpha = 0.7)
            ax.add_patch(poly)
            
            # actual calculation
            xs = t[this_idx:next_idx]
            ys = data - fit_y[-1]
            integ = np.trapz(ys, xs)
            
            
            # fig, ax2 = plt.subplots(figsize=(5,5), dpi=100)
            # ax2.plot(ts, data, '.')
            # ax2.plot(ts, fit_y)
            
            
            ### Remove point based on fit criteria ###
            # Remove if a > 0 (reductive spikes should have a < 0)
            a_idx = [idx for idx, param in enumerate(self.params) 
                     if 'a' in param]
            if popt[a_idx] > 0:
                print(f'Fitted a > 0 at {t[this_idx-delay]} s. Removing point.')
                pops.append(this_idx - delay)
                continue
            
            c_idx = [idx for idx, param in enumerate(self.params) 
                     if 'c' in param]
            if popt[c_idx] > 0:
                print(f'Fitted c > 0 at {t[this_idx-delay]} s. Removing point.')
                pops.append(this_idx - delay)
                continue 
            
            
            
            # Plot fit to main graph
            
            self.fits.append(
                            (this_idx - delay, 
                              ts+t[this_idx],
                              fit_y + baseline[this_idx:next_idx])
                             )
            
            # ax.plot(ts+t[this_idx], fit_y + baseline[this_idx:next_idx], 'gold')
            if baseline_correct:
                ax.plot(t[last_idx:next_idx], baseline[last_idx:next_idx], 'r--')
            
            fits.append(popt)
            integrals.append(integ)
            ms.append(m)
            chi_sqs.append(chi_sq)
        
        
        return fits, integrals, chi_sqs, pops, ms
    
    
    
    def analyze_file(self, baseline_correct, app_filter, 
                     delay, start_after, params, exp_func, thresh):
        '''
        Main function to call to detect and fit spikes
        '''
        
        fname = self.file.split('\\')[-1]
        print(f'\nAnalyzing {fname}')
        # Get data
        t, i, sara = self.extract_data(self.file, start_after)
        
        # Find peaks
        signals, avgFilter, stdFilter = thresholding_algo(i, lag=50, 
                                                          threshold=5,
                                                          influence = 0.6)
        # Refine peak locations
        self.idxs = self.refine_peaks(i, signals)
        
        
        # Make figure
        # if 'inline' in matplotlib.get_backend():
        #     fig, ax = plt.subplots(figsize=(5,5), dpi=300)
        # else:
        #     fig, ax = plt.subplots(figsize=(5,5), dpi=100)
        
        if app_filter:
            self.line, = self.ax.plot(t, lowpass(i,t,filter_freq))
        else:
            self.line, = self.ax.plot(t, i, color = colors[0])
        
        self.ax.set_xlabel('Time/ s')
        self.ax.set_ylabel('Current/ A')
        self.ax.set_title(fname, fontsize=9, )
        
        
        # Fit peaks. Optionally plot fits onto same ax    
        fits, areas, chi_sqs, pops, ms = self.fit_peaks(t, i, 
                                      self.idxs, sara, self.ax, 
                                      baseline_correct=baseline_correct,
                                      app_filter = app_filter, delay=delay,
                                      thresh = thresh
                                      )
        
        for pt in pops:
            # Remove spike indices that we couldn't fit
            self.idxs.pop(np.where(self.idxs==pt)[0][0])
        
        # Calculate spike integrals
        integrals, pops = self.integrate_spikes(t, i, self.idxs, avgFilter)
        
        for pt in pops:
            # Remove spike indices with negative integrals
            self.idxs.pop(np.where(self.idxs==pt)[0][0])
        
        # Get time points
        points = t[self.idxs]
        
        self.draw_fits()
        
        # Transpose data for saving
        fits = np.array(fits).T
        
        self.save_fits  = fits
        self.areas      = areas
        self.integrals  = integrals
        self.chi_sqs    = chi_sqs
        self.points     = points
        self.ms         = ms
        
        return fits, areas, integrals, chi_sqs, points, ms



class Index:
        
    def __init__(self, folder, fig, ax, params, exp_func):
        
        self.ind   = 0
        self.folder= folder
        self.files = [os.path.join(folder, file)
                      for file in os.listdir(folder)
                      if file.endswith('.txt')]
        self.inds  = [i for i in range(len(self.files))]
        self.file  = self.files[0]
        
        self.fig   = fig
        self.ax    = ax
        
        self.Picker = {}
        self.data   = {}
        
        self.params = params
        self.exp_func = exp_func
        
                
        self.Picker[self.ind] = InteractivePicker(self.file, self.fig, self.ax,
                                                  self.params, self.exp_func)
        self.Picker[self.ind].analyze_file(correct_baselines, apply_filter, 
                                           delay, start_after, self.params,
                                           self.exp_func, thresh)
        
        self.load_cids(self.Picker[self.ind])
        
        plt.show()
    
    
    def load_cids(self, Picker):
        if hasattr(self, 'cid'):
            for cid in [self.cid, self.cid2, self.cid3]:
                self.fig.canvas.mpl_disconnect(cid)
        self.cid  = self.fig.canvas.mpl_connect('button_press_event', Picker)
        self.cid2 = self.fig.canvas.mpl_connect('key_press_event', Picker.keypress_event_handler)
        self.cid3 = self.fig.canvas.mpl_connect('scroll_event', Picker.scroll_event_handler)
      
        
    def _next(self, event):
        
        self.store_data()
        self._del_picker(self.Picker[self.ind])
        
        self.ind += 1
        if self.ind > self.inds[-1]:
            self.ind = self.inds[0]
        
        self.file = self.files[self.ind]
        
        self.ax.clear()
        
        self.Picker[self.ind] = InteractivePicker(self.file, self.fig, self.ax,
                                                  self.params, self.exp_func)
        self.Picker[self.ind].analyze_file(correct_baselines, apply_filter, 
                                           delay, start_after, self.params,
                                           self.exp_func, thresh)
        
        self.load_cids(self.Picker[self.ind])
        
        plt.show()
        
    
    def _prev(self, event):
        
        self.store_data()
        self._del_picker(self.Picker[self.ind])
        
        self.ind -= 1
        if self.ind < self.inds[0]:
            self.ind = self.inds[-1]
            
        self.file = self.files[self.ind]
        
        self.ax.clear()
        
        self.Picker[self.ind] = InteractivePicker(self.file, self.fig, self.ax,
                                                  self.params, self.exp_func)
        self.Picker[self.ind].analyze_file(correct_baselines, apply_filter, 
                                           delay, start_after, self.params,
                                           self.exp_func, thresh)
        
        self.load_cids(self.Picker[self.ind])
        
        plt.show()
    
        
    def _del_picker(self, picker):
        for attr, _ in list(picker.__dict__.items()):
            delattr(picker, attr)
    
    
    def reset(self, event):
        self.Picker[self.ind]._reset(event)
     
     
    def store_data(self):
        
        d = {}
        d['Index']          = []
        d['time/ s']        = []
        d.update({param: [] for param in self.params})
        d['Catalytic area/ C'] = []
        d['baseline slope'] = []
        d['integral/ C']    = []
        d['Reduced Chi^2']  = []
        d['File']           = []
        
        p = self.Picker[self.ind]
        
                
        for i in range(len(p.integrals)):
            d['Index'].append(i+1)
        for i in range(len(self.params)):
            for x in p.save_fits[i]:
                d[params[i]].append(x)
        for x in p.areas:
            d['Catalytic area/ C'].append(x)
        for x in p.ms:
            d['baseline slope'].append(x)
        for x in p.integrals:
            d['integral/ C'].append(x)
        for x in p.chi_sqs:
            d['Reduced Chi^2'].append(x)
        for x in p.points:
            d['time/ s'].append(x)
        
        d['File'].append(p.file)
        d['File'].append(f'Fit: {fit}')
        d['File'].append(f'Baseline correct: {correct_baselines}')
        d['File'].append(f'Delay: {delay} pts')
        
        for key in d:
            while len(d[key]) < max([len(d[key]) for key in d]):
                d[key].append(' ')
    
        for key in d:
            d[key].append(' ')
        
        
        self.data[self.ind] = d
                
                
     
    def save(self, event):
        
        self.store_data()
        
        dout  = {}
        dicts = [self.data[key] for key in self.data.keys()]
        for key in dicts[0].keys():
            dout[key] = list(d[key] for d in dicts)
            # Magic one liner
            # for sublist in list: for item in sublist: yield item
            dout[key] = [item for sublist in dout[key] for item in sublist]
        
        
        df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dout.items() ]))

    
        writer = pd.ExcelWriter(self.folder + '/output.xlsx', engine='xlsxwriter')
        df.to_excel(writer, index=False, header=True, startcol=0)
        writer.save()
        print(f'Saved as {self.folder}/output.xlsx')
        
        return


def shift_axes(ax, direction):
    # Pan axes if arrow keys are pressed
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xdelta = 0.2*(xlim[1] - xlim[0])
    ydelta = 0.2*(ylim[1] - ylim[0])
    if direction == 'right':
        xlim = [
            xlim[0] + xdelta,
            xlim[1] + xdelta,
            ]
    if direction == 'left':
        xlim = [
            xlim[0] - xdelta,
            xlim[1] - xdelta,
            ]
    if direction == 'up':
        ylim = [
            ylim[0] + ydelta,
            ylim[1] + ydelta,
            ]
    if direction == 'down':
        ylim = [
            ylim[0] - ydelta,
            ylim[1] - ydelta,
            ]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return


def zoom_axes(ax, direction, center):
    # Zoom axes by mouse scroll wheel
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    
    # Shift (almost) to new center
    x, y = center
    shift_x = x - (xlim[0] + xlim[1])/2 
    shift_y = y - (ylim[0] + ylim[1])/2 
    
    if direction == 'down':
        shift_x = -shift_x
        shift_y = -shift_y
    
    xlim = [
        xlim[0] + 0.7*shift_x,
        xlim[1] + 0.7*shift_x,
        ]  
    ylim = [
        ylim[0] + 0.7*shift_y,
        ylim[1] + 0.7*shift_y,
        ]
        
    # Zoom in or out
    xdelta = 0.2*(xlim[1] - xlim[0])
    ydelta = 0.2*(ylim[1] - ylim[0])
    
    if direction == 'down':
        xdelta = -xdelta
        ydelta = -ydelta
    
    xlim = [
        xlim[0] + xdelta,
        xlim[1] - xdelta,
        ]
    ylim = [
        ylim[0] + ydelta,
        ylim[1] - ydelta,
        ]
        
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return        
    


if __name__ == '__main__':

    fig, ax = plt.subplots(figsize=(5,6), dpi=100)    
    plt.subplots_adjust(bottom=0.3)
    
    index = Index(data_folder, fig, ax, params, exp_func)
    
    
    axcalc = plt.axes([0.4, 0.1, 0.25, 0.05])
    bcalc = Button(axcalc, 'Reset')
    bcalc.on_clicked(index.reset)   
    
    axprev = plt.axes([0.15, 0.1, 0.25, 0.05])
    bnext = Button(axprev, 'Previous')
    bnext.on_clicked(index._next) 
    
    axnext = plt.axes([0.65, 0.1, 0.25, 0.05])
    bprev = Button(axnext, 'Next')
    bprev.on_clicked(index._prev) 
    
    axsave = plt.axes([0.4, 0.025, 0.25, 0.05])
    bsave = Button(axsave, 'Save')
    bsave.on_clicked(index.save)
