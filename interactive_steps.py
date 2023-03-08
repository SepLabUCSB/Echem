import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from matplotlib.widgets import Button
from matplotlib.widgets import Slider
from matplotlib.widgets import CheckButtons
from matplotlib.widgets import TextBox
import os
# matplotlib.use('qt5agg')
plt.ion() # Interactive mode
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

# folder = r'C:\Users\Eric Liu\Desktop\LiorLab\InsulatingProject\cubese12'
# folder = r'C:\Users\BRoehrich\Desktop\data'
folder = r'C:\Users\BRoehrich\Desktop\Julia_data'

"""
**Type %matplotlib in console before running file**

Python script for interactive analysis of single entity echem step data.

Analyzes all files in specified folder.

The program will automatically try to find steps (with >2% 

relative change). Click to select or deselect step points.

Click "Export" when finished to save a list of all step sizes

to a folder of your choice.

NOTE: the script assumes current data is given in mA. Edit

the get_data function to change this.

"""


# PARAMETERS TO ADJUST
POTENTIOSTAT = 'Biologic'
# POTENTIOSTAT = 'HEKA'

threshold  = 1     # Autopick steps with 1st derivative n stdevs above average
cutoff     = 10    # Remove first n seconds




filtered   = False    # Apply Savinsky-Golay filter to plotted current (raw data unfiltered)
filter_win = 101      # S-G filter window (int)
filter_ord = 1        # S-G filter order (int)

fit_type   = 'savgol' # Fit to do between steps
linear_fit = True     # Do linear fit between points. If False, uses median
                      #    value between points instead. True accounts for 
                      #    non-zero baseline. Either way, the values are the 
                      #    raw, non-filtered data!








def get_data(file, potentiostat):
    if file.endswith('.asc'): 
        potentiostat = 'HEKA'
        
    if potentiostat == 'Biologic':
        df = pd.read_fwf(file, skiprows=1, headers=0,
                              names=('t', 'i'))
        df['i'] = df['i']/1000 #convert mA -> A
        df = df[df['t'] > cutoff]
        
        # df = df.groupby(np.arange(len(df))//10).mean() #compress to 100 Hz
    
        return df
    
    
    elif potentiostat == 'HEKA':        
        # Use StringIO for faster file reading, and handling multiple
        # PATCHMASTER Series in a single file
        from io import StringIO
        
        s = StringIO()
        
        def isfloat(x):
            try:
                float(x)
                return True
            except:
                return False
        
        with open(file) as f:
            for line in f:
                if isfloat(line.split(',')[0]):
                    #skip rows which don't start with the index number (1, 2, ...)
                    s.write(line)
    
        s.seek(0) #return to top of StringIO
        
        
        df = pd.read_csv(s, names=(
            "index", "t", "i", "time2", "v"), 
            engine='c', dtype=float)
        df = df.drop(['time2', 'index', 'v'], axis=1)
        df = df[df['t'] > cutoff]
        
        if df.empty:
            df = get_data(file, 'Biologic')
        
        return df


def isEven(x):
    return not bool(x%2)


class StepPicker:
    
    global threshold, linear_fit, filtered, filter_win, filter_ord
    
    """Based on
    https://scipy-cookbook.readthedocs.io/items/Matplotlib_Interactive_Plotting.html
    """

    def __init__(self, xdata, ydata, plot_ydata, points=[], ax=None):
        self.xdata = np.array(xdata)            # Time
        self.ydata = np.array(ydata)            # Raw current data
        self.plot_ydata = np.array(plot_ydata)  # (possibly) filtered current
        self.criticalPoints = dict()            # Step points (x,y):1
        self.steps = []
        self.steptimes = []
        
        self.frequency = 1/np.mean(np.diff(xdata))
        
        self.min_spacing = max(5, int(0.2 * self.frequency)) # at least 5 points, or 0.2 s
        # print(f'freq: {self.frequency:0.2f}, spacing: {self.min_spacing}')
        
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax
        
        # self.line = average (stairstep) line
        # self.drawnpoints = step location markers
        self.line, = self.ax.plot([self.xdata[0]], [self.ydata[0]], '.-') 
        self.drawnpoints = self.ax.add_line(
            matplotlib.lines.Line2D([],[],linewidth=0, 
                                    marker='d', c='r', zorder=100))
        
          
        foundPoints = self.detect_steps()
        for (x,y) in foundPoints:
            self.drawPoints(ax, x, y)
        
        self.ax.figure.canvas.draw_idle()


    def _update(self, xdata=None, ydata=None, 
                pts=None, ax=None, avg=None):
        '''
        Call to redraw avg line and critical point
        markers. Can pass new x,y data, add critical points, etc.

        '''
        if xdata is not None:
            self.xdata = xdata
            
        if ydata is not None:
            self.ydata = ydata
            
        if ax is not None:
            self.ax = ax
        
        if avg is not None:
            self.avg = avg
        


        # Redraw avg line and critical points
        try:
            self.line, = self.ax.plot([self.xdata[0]], 
                                      [self.ydata[0]], color=colors[0]) 
            self.line.set_data(self.xdata, self.avg)
            self.ax.figure.canvas.draw_idle()
        except:
            # Exception if self.avg hasn't been generated yet
            pass


        self.drawnpoints = self.ax.add_line(matplotlib.lines.Line2D(
            [],[],linewidth=0, marker='d', c='r', zorder=100))
        
        xlist = []
        ylist = []
        for (x,y) in self.criticalPoints:
            xlist.append(x)
            ylist.append(y)
        
        self.drawnpoints.set_data(xlist, ylist)
            
        self.ax.figure.canvas.draw_idle()

            

    
    def __call__(self, event):
        '''
        Called on mouse click in graph axes
        
        Prioritize removing marked point near cursor
        (+- xtol seconds)
        
        Otherwise add a new point at (x, y(x))
        '''
        
        # Set xtol based on current axis limits 
        upper = self.ax.get_xlim()[1]
        lower = self.ax.get_xlim()[0]
        diff = upper-lower
        xtol = 0.01*diff
        
        deleted = False
        if event.inaxes and fig.canvas.manager.toolbar.mode == "":
            # print(xtol)
            clickX = event.xdata
            if (self.ax is None) or (self.ax is event.inaxes):
                # Prioritizing removing marked point
                for (x,y) in list(self.criticalPoints):
                    if (clickX-xtol < x < clickX+xtol):
                        self.drawPoints(event.inaxes, x, y)
                        deleted = True
                    else:
                        continue
                # Otherwise, add a new point
                if not deleted:
                    i = np.abs(self.xdata-clickX).argmin()
                    this_x = self.xdata[i]
                    this_y = self.plot_ydata[i]
                    self.drawPoints(event.inaxes, this_x, this_y)
    
    
    def keypress_event_handler(self, event):
        key = event.key
        if key in ['left', 'right', 'up', 'down']:
            shift_axes(self.ax, key)
        self.ax.figure.canvas.draw_idle()
    
    
    def scroll_event_handler(self, event):
        key = event.button
        x, y = event.xdata, event.ydata
        if key in ['down', 'up']:
            zoom_axes(self.ax, key, (x,y))
        self.ax.figure.canvas.draw_idle() 

    
    def drawPoints(self, ax, x, y):
        """
        Draw or remove the point on the plot and from self.criticalPoints
        """
        if (x, y) in self.criticalPoints:
            self.criticalPoints.pop((x,y), None)
            
            xlist = []
            ylist = []
            for (x,y) in self.criticalPoints:
                xlist.append(x)
                ylist.append(y)
                
            self.drawnpoints.set_data(xlist, ylist)
            self.ax.figure.canvas.draw_idle()
            # print(len(self.criticalPoints))

            
        else:
            # Draw new point, add to criticalPoints dict

            self.criticalPoints[(x, y)] = 1
            xlist = []
            ylist = []
            for (x,y) in self.criticalPoints:
                xlist.append(x)
                ylist.append(y)
            
            self.drawnpoints.set_data(xlist, ylist)
                
            self.ax.figure.canvas.draw_idle()
            # print(len(self.criticalPoints))
    
    
    def remove_hanging_points(self, xmin, xmax):
        '''
        Removes points which are outside of the currently selected 
        time boundaries xmin: xmax
        '''
        for (x, y) in self.criticalPoints.copy():
            if (x < xmin) or (x > xmax):
                # Remove point
                self.drawPoints(self.ax, x, y)
        
        
    
    def detect_steps(self, thresh=threshold):
        '''
        Initial step detection algorithm
        
        Apply savitzky-golay filter, calculate 1st derivative,
        pick steps based on spikes in 1st derivative
        '''
        
        window = 3*self.min_spacing
        if isEven(window):
            window -= 1
        dy = savgol_filter(self.ydata, window, 1, deriv=1, mode='nearest')
        self.dy = dy
        dy = abs(dy)
        
        avg = np.average(dy)
        std = thresh*np.std(dy)
        
        
        # fig, ax = plt.subplots(figsize=(5,5), dpi=100)
        # ax.plot(dy)
        # ax.plot(np.ones(len(dy))*avg)
        # ax.plot(np.ones(len(dy))*(avg+std), 'orange')
        # ax.plot(np.ones(len(dy))*(avg-std), 'orange')
        
        # idxs = np.where(abs(dy) > avg + std)
        
        idxs = []
        in_spike = False
        _out_counter = 0
        
        # Identify spikes by changes in velocity
        for i in range(len(dy)):
            
            if (abs(dy[i]) >= avg + std and not in_spike):
                # First point above std
                idxs.append(i)
                in_spike = True
            
            elif (abs(dy[i]) >= avg + std and in_spike):
                # subsequent point in same spike
                continue
            
            elif (abs(dy[i]) < avg + std and in_spike):
                # second crossing point
                _out_counter += 1
                if _out_counter > self.min_spacing:
                    in_spike = False
                    _out_counter = 0
            
        # Verify steps have a reasonable change in current
        # for i in idxs[:]:
        #     window = 50
            
        #     val_before = np.average(data[i-window:i-5])
        #     val_after  = np.average(data[i+5:i+window])
            
        #     delta = val_after - val_before
            
        #     if abs(delta/val_before) < 0.0015:
        #         idxs.remove(i)
        
        points = [(self.xdata[i], self.plot_ydata[i]) for i in idxs]
        
        return points
        
    
    
    def calculate_steps(self, button):
        '''
        Refines step locations, then draws fit line and saves steps
        '''
        delta = np.zeros(len(self.xdata))
        self.avg = np.zeros(len(self.xdata))
        
        
        # Refine step locations
        for i in range(len(self.xdata)-1):
            delta[i] = abs((self.plot_ydata[i+1]-self.plot_ydata[i])/self.plot_ydata[i])
        
        # delta = self.dy
        
        
        for (x,y) in sorted(list(self.criticalPoints), key=lambda x:x[0]):
            
            # find largest local step, search +- m points
            m = self.min_spacing
            xi = np.where(self.xdata == x)[0][0] #convert to index
            
                
            n = np.where(delta == max(delta[xi-m:xi+m]))[0][0]
            if not n == xi:
                # Remove old point
                # print(f'Removing {x:0.3f}. Found better point {self.xdata[n]:0.3f}')
                self.drawPoints(self.ax, x, y)
                # Add new point (xdata[n], ydata[n])
                if (self.xdata[n], self.plot_ydata[n]) not in self.criticalPoints:
                    self.drawPoints(self.ax, self.xdata[n], 
                                    self.plot_ydata[n])
        
        
        # Smooth data between steps to measure step size
        # more accurately
        indices = [np.where(self.xdata == x)[0][0] 
                   for (x,y) in list(self.criticalPoints)]
    
        indices.insert(0,0)
        indices.append(len(self.xdata-1))
        indices.sort()
        noises = []
        
        for i in range(len(indices)-1):
            index = indices[i]
            next_index = indices[i+1]
            
            pad = (next_index - index)//10
            if pad*self.frequency > 0.5:
                pad = int(0.5/self.frequency)
            if pad < 2: pad = 2
            
            this_ydata = self.ydata[index+pad: next_index-pad]
            window = len(this_ydata)//5 
            if window == 1:
                window = 3
            if isEven(window):
                window += 1
            if window > 100:
                window = 99
            
            if window == 1:
                filtered = this_ydata
            else:
                filtered = savgol_filter(this_ydata, window,
                                     polyorder=1)
            
            # Calculate noise as data - filtered data
            this_noise = np.sqrt(
                np.mean( (this_ydata-filtered)**2 )
                )
            noises.append(this_noise)
            filtered = np.pad(filtered, pad, mode='constant',
                              constant_values=(this_ydata[0],
                                               this_ydata[-1]))
            self.avg[index:next_index] = filtered
                        
            # if linear_fit:
            #     # Fit line in between steps to account for sloped baseline
            #     m, b = np.polyfit(np.arange(index+pad,next_index-pad), 
            #                           this_ydata, 1)
            
            #     for i in range(index, next_index):
            #         self.avg[i] = m*i + b
                
            #     this_noise = np.sqrt(np.mean((this_ydata - self.avg[index+pad:next_index-pad])**2))
            #     noises.append(this_noise)
                    
            
            # else:
            #     medianvalue = np.median(this_ydata)
            #     for i in range(index, next_index):
            #         self.avg[i] = medianvalue
                
            
            
        
            
                
            if index != 0:
                self.avg[index] = self.avg[index-1]
             
        # Draw result on graph and save step sizes 
        self.draw_average()
        self.get_steps()
        self.noise = np.average(noises[2:len(noises)-2])
        print(f'Average RMS noise: {self.noise: .1E} A')
    
        
    
            
    def draw_average(self):
        '''
        Redraw smoothed step line
        '''
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        self.line.remove()
        self.line, = self.ax.plot([self.xdata[0]], 
                                  [self.ydata[0]], 
                                  color=colors[0]) 
        
        self.line.set_data(self.xdata, self.avg)
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        self.line.figure.canvas.draw()
    
    
    
    
    def get_steps(self):
        '''
        Save current step sizes to self.steps and self.abs_steps
        '''
        self.steptimes = []
        self.steps = []
        self.abs_steps = []
        for (x,y) in sorted(list(self.criticalPoints), 
                            key=lambda x:x[0]):
            i = np.where(self.xdata==x)[0][0]
            step = (self.avg[i+1]-self.avg[i])/self.avg[i]
            abs_step = self.avg[i+1]-self.avg[i]
            self.steptimes.append(x)
            self.steps.append(step)
            self.abs_steps.append(abs_step)




class Index:
    '''
    Class for cycling through multiple graphs
    '''
    
    global folder, files, checkbox, pointsbox, POTENTIOSTAT, filtered, filter_ord, filter_win
    
    
    def __init__(self):
        

        self.ind = 0
        self.sp = dict() # sp[i] is StepPicker obj
        self.files = files
        
        self.slider  = dict()
        self.slider2 = dict()
        self.xs      = dict()
        self.ys      = dict()
        self.plot_ys = dict()
        
        i = self.ind % len(files)
        self.i = i
        self.slider[i] =  0 #initialize slider index
        
        
        df = get_data(files[i], POTENTIOSTAT)
        
        # Initialize plot
        name = files[i].split('/')[-1][:-4]
        ax.set_title('%s: %s'%((i+1), name), pad=15)
        
        self.xs[i] = df['t'].to_numpy()
        self.ys[i] = df['i'].to_numpy()
        
        if filtered:
            from scipy.signal import savgol_filter
            self.plot_ys[i] = savgol_filter(self.ys[i], 
                                       filter_win, filter_ord)
        else:
            self.plot_ys[i] = self.ys[i]
        
        self.line, = ax.plot(self.xs[i], self.plot_ys[i], 'k-')
        self.sp[i] = StepPicker(self.xs[i], self.ys[i], 
                                self.plot_ys[i], ax=ax)
        
        self.cid = fig.canvas.mpl_connect('button_press_event', 
                                          self.sp[i])
        self.cid2 = fig.canvas.mpl_connect('key_press_event',
                                           self.sp[i].keypress_event_handler)
        self.cid3 = fig.canvas.mpl_connect('scroll_event',
                                           self.sp[i].scroll_event_handler)
        
        self.slider2[i] = len(self.xs[i]) - 1
        
        plt.show()
    
        
    
    def next(self, event):
        '''
        Load next file
        '''
        ax.clear()
        fig.canvas.mpl_disconnect(self.cid)
        fig.canvas.mpl_disconnect(self.cid2)
        self.ind += 1
        i = self.ind % len(files)  
        self.i = i
        
        name = files[i].split('/')[-1][:-4]
        ax.set_title('%s: %s'%((i+1), name), pad=15)
        
        if i not in self.sp:
            # Create new
            df = get_data(files[i], POTENTIOSTAT)
            
            self.xs[i] = df['t'].to_numpy()
            self.ys[i] = df['i'].to_numpy()
        
            if filtered:
                from scipy.signal import savgol_filter
                self.plot_ys[i] = savgol_filter(self.ys[i], 
                                           filter_win, filter_ord)
            else:
                self.plot_ys[i] = self.ys[i]
            
            slider.set_val(self.xs[i][0])
            slider2.set_val(self.xs[i][-1])
            self.line, = ax.plot(self.xs[i], self.plot_ys[i], 'k-')
            self.sp[i] = StepPicker(self.xs[i], self.ys[i], 
                                    self.plot_ys[i], ax=ax)
            
        else:
            # Reinitialize
            
            ind  = self.slider[i]
            ind2 = self.slider2[i]
            slider.set_val(self.xs[i][ind])
            slider2.set_val(self.xs[i][ind2])
            
            self.line, = ax.plot(self.xs[i][ind:ind2], 
                                 self.plot_ys[i][ind:ind2], 'k-')
            
            self.sp[i]._update()
            
        self.cid = fig.canvas.mpl_connect('button_press_event', self.sp[i])
        self.cid2 = fig.canvas.mpl_connect('key_press_event',
                                           self.sp[i].keypress_event_handler)
        
        plt.show()
        
    
    def prev(self, event):
        '''
        Load previous file
        '''
        ax.clear()
        fig.canvas.mpl_disconnect(self.cid)
        fig.canvas.mpl_disconnect(self.cid2)
        self.ind -= 1
        i = self.ind % len(files)
        self.i = i
        
        name = files[i].split('/')[-1][:-4]
        ax.set_title('%s: %s'%((i+1), name), pad=15)
        
        if i not in self.sp:
            # Create new
            df = get_data(files[i], POTENTIOSTAT)
            
            self.xs[i] = df['t'].to_numpy()
            self.ys[i] = df['i'].to_numpy()
        
            if filtered:
                from scipy.signal import savgol_filter
                self.plot_ys[i] = savgol_filter(self.ys[i], 
                                           filter_win, filter_ord)
            else:
                self.plot_ys[i] = self.ys[i]
            
            slider.set_val(self.xs[i][0])
            slider2.set_val(self.xs[i][-1])
            self.line, = ax.plot(self.xs[i], self.plot_ys[i], 'k-')
 
            self.sp[i] = StepPicker(self.xs[i], self.ys[i], 
                                    self.plot_ys[i], ax=ax)
            
        else:
            # Reinitialize
            
            ind  = self.slider[i]
            ind2 = self.slider2[i]
            slider.set_val(self.xs[i][ind])
            slider2.set_val(self.xs[i][ind2])
            
            self.line, = ax.plot(self.xs[i][ind:ind2], 
                                 self.plot_ys[i][ind:ind2], 'k-')
            
            self.sp[i]._update()
            
        self.cid = fig.canvas.mpl_connect('button_press_event', self.sp[i])
        self.cid2 = fig.canvas.mpl_connect('key_press_event',
                                           self.sp[i].keypress_event_handler)
    
        
        plt.show()        
    
        
    def recalc(self, event):
        '''
        Call calculate_steps on current graph instance
        
        Accounts for slider removing first data points
        '''
        # i = self.ind % len(files)
        i = self.i
        ind  = self.slider[i]
        ind2 = self.slider2[i]
        
        try:
            self.sp[i].xdata = self.xs[i][ind:ind2]
            self.sp[i].ydata = self.ys[i][ind:ind2]
            self.sp[i].plot_ydata = self.plot_ys[i][ind:ind2]
            
            # Remove points outside slider selected window
            xmin, xmax = self.xs[i][ind], self.xs[i][ind2]
            self.sp[i].remove_hanging_points(xmin, xmax)

            self.sp[i].calculate_steps(event)
        
        except IndexError:
            print('ERROR: remove hanging points')
    
    
    def hist(self, event):
        '''
        Display current histogram. Absolute value
        
        checkbox determines whether relative or absolute
        
        step sizes are used
        '''
        step_list = []
        
        abs_deltaI = checkbox.get_status()
        
        if abs_deltaI[0] is False:
            for i in range(len(files)):
                try:
                    steps = self.sp[i].steps
                    n = 0
                    for step in steps:
                        if n < int(pointsbox.text):
                            step_list.append(step)
                            n += 1
                except KeyError:
                    print('No steps saved in File %s' %(i+1))
            
            xlabel = '$\Delta$$I/I_{ss}$'
            
        if abs_deltaI[0] is True:
            for i in range(len(files)):
                try:
                    steps = self.sp[i].abs_steps
                    n = 0
                    for step in steps:
                        if n < int(pointsbox.text):
                            step_list.append(step)
                            n += 1
                except KeyError:
                    print('No steps saved in File %s' %(i+1))
           
            xlabel = '$\Delta I/$ $A$'
            
            
        if step_list:
            plt.figure(figsize=(6,6), dpi=100)
            # bins = np.arange(0, 1, 0.005)
            plt.hist(step_list, rwidth=0.8)
            # plt.xlim(-0.005, 1.1*max(step_list))
            plt.xlabel(xlabel)
            plt.ylabel('Count')
            plt.text(1, 1.02, 'N = %s' %len(step_list),
                     transform=ax.transAxes,
                     horizontalalignment='right')
            plt.show()
        
         
    
    
    def save(self, event):
        '''
        Export list of relative and absolute step
        
        sizes as xlsx document
        '''
        file_index = []
        file_list = []
        time_list = []
        step_list = []
        abs_step_list = []
        noise_list = []
        
            
        for i in range(len(files)):
            # Create list of steps and abs_steps to export
            try:
                file = self.files[i]
                times = self.sp[i].steptimes
                steps = self.sp[i].steps
                abs_steps = self.sp[i].abs_steps
                noise = self.sp[i].noise
                
                n = 0
                m = 0

                for time, step in zip(times, steps):
                    if n < int(pointsbox.text):
                        # Only keep first n steps
                        file_index.append(i+1)
                        file_list.append(file)
                        time_list.append(time)
                        step_list.append(step)
                        noise_list.append(noise)
                        n += 1
                for abs_step in abs_steps:
                    if m < int(pointsbox.text):
                        abs_step_list.append(abs_step)
                        m += 1
                
                file_index.append('')
                file_list.append('')
                time_list.append('')
                step_list.append('')
                abs_step_list.append('')
                noise_list.append('')
                
                        
            except KeyError:
                print('No steps saved in File %s' %(i+1))
        
        
        if step_list:
            steps_file = folder + '/steps.xlsx'
            
            writer = pd.ExcelWriter(steps_file, engine = 'xlsxwriter')
            out = pd.DataFrame(
                {'File #': file_index,
                'File': file_list,
                'Time/s': time_list,
                'dI/Iss': step_list,
                'delta I (A)': abs_step_list,
                'RMS noise (A)': noise_list
                }
                )
            out.to_excel(writer, index=False, header=True, startcol=0)
            writer.save()
            print('Saved as %s' %steps_file)
        else:
            print('No steps!')
    
    
    
    def slider_changed(self, val):
        '''
        Redraw raw data after adding/ removing xdata using slider
        '''        
        
        i = self.i
        
        def find_nearest(array,value):
            idx = np.searchsorted(array, value, side="left")
            if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
                return array[idx-1]
            else:
                return array[idx]
        
        try:
            ind = np.where(self.xs[i] == find_nearest(self.xs[i], val))[0][0]
            self.slider[i]  = ind
            ind2 = self.slider2[i]
            
            self.line.set_xdata(self.xs[i][ind:ind2])
            self.line.set_ydata(self.plot_ys[i][ind:ind2])
            ax.set_xlim(0.95*self.xs[i][ind],
                        1.05*self.xs[i][ind2])
            ax.figure.canvas.draw_idle()

        except KeyError:
            pass
    
    
    def slider2_changed(self, val):
        '''
        Redraw raw data after adding/ removing xdata using slider
        '''        
        
        i = self.i
        
        def find_nearest(array,value):
            idx = np.searchsorted(array, value, side="left")
            if idx > 0 and (idx == len(array) or abs(value - array[idx-1]) < abs(value - array[idx])):
                return array[idx-1]
            else:
                return array[idx]
        
        try:
            ind2 = np.where(self.xs[i] == find_nearest(self.xs[i], val))[0][0]
            self.slider2[i]  = ind2
            ind = self.slider[i]
            
            self.line.set_xdata(self.xs[i][ind:ind2])
            self.line.set_ydata(self.plot_ys[i][ind:ind2])
            ax.set_xlim(0.95*self.xs[i][ind],
                        1.05*self.xs[i][ind2])
            ax.figure.canvas.draw_idle()

        except KeyError:
            pass
        
        
  

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
        xdelta = - xdelta
        ydelta = - ydelta
    
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
        




# Unbind left and right arrow keys so we can use them to pan the graph
try:
    plt.rcParams['keymap.back'].remove('left')
    plt.rcParams['keymap.forward'].remove('right')
except: # Already removed
    pass


os.chdir(folder)
files = []
for file in os.listdir(folder):
    if (file.endswith('.txt')
        or file.endswith('.asc')):
        files.append(os.path.join(folder,file))


fig, ax = plt.subplots(figsize=(5,6), dpi=100)
plt.subplots_adjust(bottom=0.3)


callback = Index()

# Recalculate step sizes
axcalc = plt.axes([0.5, 0.075, 0.25, 0.05])
bcalc = Button(axcalc, 'Recalculate')
bcalc.on_clicked(callback.recalc)

# Next file
axnext = plt.axes([0.8, 0.075, 0.1, 0.05])
bnext = Button(axnext, 'Next')
bnext.on_clicked(callback.next)

# Previous file
axprev = plt.axes([0.35, 0.075, 0.1, 0.05])
bprev = Button(axprev, 'Prev')
bprev.on_clicked(callback.prev)

# Save as xlsx
axexport = plt.axes([0.1, 0.075, 0.2, 0.05])
bexport = Button(axexport, 'Export')
bexport.on_clicked(callback.save)

# Plot histogram
axplotbutton = plt.axes([0.1, 0.005, 0.3, 0.05])
histbutton = Button(axplotbutton, 'Plot histogram')
histbutton.on_clicked(callback.hist)

# Starting point slider
axslider = plt.axes([0.1, 0.2, 0.8, 0.025])
slider = Slider(axslider, '', callback.xs[0][0], 
                callback.xs[0][-1], valinit=callback.xs[0][0])
slider.on_changed(callback.slider_changed)

# Ending point slider
axslider2 = plt.axes([0.1, 0.15, 0.8, 0.025])
slider2 = Slider(axslider2, '', callback.xs[0][0], 
                callback.xs[0][-1], valinit=callback.xs[0][-1])
slider2.on_changed(callback.slider2_changed)

# Absolute delta I checkbox
axcheckbox = plt.axes([0.4, 0.005, 0.35, 0.05])
checkbox = CheckButtons(axcheckbox, ['Absolute $\Delta$I'])

# Select first n points box
axpointsbox = plt.axes([0.8, 0.005, 0.1, 0.05])
pointsbox = TextBox(axpointsbox, '', initial='10')



