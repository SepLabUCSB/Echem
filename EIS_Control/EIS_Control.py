import tkinter as tk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from cycler import cycler
import os
import sys
import time
import pyvisa
import rigol_control
import siglent_control
import create_waveform

this_dir = rigol_control.__file__[:-16]


# Find dirs with waveforms
csv_dir = os.path.join(this_dir, 'csv')
rigol_waves = os.path.join(this_dir, 'rigol_waves')


# Get matplotlib style sheet and color cycle
plt.style.use(os.path.join(this_dir, 'scientific.mplstyle'))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']



class MainWindow:
    
    global csv_dir, rigol_waves
    
    def __init__(self, root):
        self.root = root
        root.title("Arb/ Scope control")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    

        self.rm = pyvisa.ResourceManager()
        self.ft = None
            
        
        # Initialize frames and canvas
        self.frame = tk.Frame(self.root)
        self.frame.grid(row=0, column = 0)
        
        self.frame2 = tk.Frame(self.root)
        self.frame2.grid(row=0, column = 1)
        
       
        self.fig = plt.Figure(figsize=(4,4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        # self.ax.plot([],[])
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=1, column=0)
        
        
        
        #################################################
        ### FRAME 1: Instrument selection and control ###
        #################################################
        
        # Waveform selection dropdown menu
        text = tk.Label(self.frame, text='Waveform:')
        text.grid(row=2, column=1)
        self.file_list = [file for file in os.listdir(rigol_waves) 
                          if file.endswith('freqs.csv')
                          if file.startswith('Rigol')]
        
        self.waveform = tk.StringVar(self.frame)
        self.waveform.set(self.file_list[0])
        self.waveform_selector = tk.OptionMenu(self.frame, self.waveform, 
                                               *self.file_list, command=self.show_waveform)
        self.waveform_selector.grid(row=2, column=2)

                
        
        # VISA selection menu: Arb
        text = tk.Label(self.frame, text='Arb:')
        text.grid(row=3, column=1)
        self.arb = tk.StringVar(self.frame)
        self.arb.set('USB0::0x1AB1::0x0643::DG8A232302748::INSTR')
        self.arb_selector = tk.OptionMenu(self.frame, self.arb, 
                                               *self.rm.list_resources())
        self.arb_selector.grid(row=3, column=2)
        self.apply_waveform_button = tk.Button(self.frame, text='Apply Wave', 
                                               command=self.apply_waveform)
        self.apply_waveform_button.grid(row=3, column=3)
        
        
        
        # VISA selection menus: Scope
        text = tk.Label(self.frame, text='Scope:')
        text.grid(row=4, column=1)
        self.scope = tk.StringVar(self.frame)
        self.scope.set('USB0::0xF4ED::0xEE3A::SDS1EDEX5R5381::INSTR')
        self.scope_selector = tk.OptionMenu(self.frame, self.scope, 
                                               *self.rm.list_resources())
        self.scope_selector.grid(row=4, column=2)
        self.record_signals_button = tk.Button(self.frame, text='Record Signals', 
                                               command=self.record_signals)
        self.record_signals_button.grid(row=4, column=3)
        
        
        
        ########################################
        ### FRAME 2: User-adjustable options ###
        ########################################
        

        # Applied waveform amplitude
        text = tk.Label(self.frame2, text='Waveform Vpp (mV):')
        text.grid(row=0, column = 0)
        self.waveform_vpp = tk.Text(self.frame2, height=1, width=7)
        self.waveform_vpp.insert('1.0', '20')
        self.waveform_vpp.grid(row=0, column=1)
        
        
        
        # Recording duration
        text = tk.Label(self.frame2, text='Recording time (s):')
        text.grid(row=1, column = 0)
        self.recording_time = tk.Text(self.frame2, height=1, width=7)
        self.recording_time.insert('1.0', '10')
        self.recording_time.grid(row=1, column=1)
        
        
        
        # Autolab current range
        text = tk.Label(self.frame2, text='Current range:')
        text.grid(row=2, column = 0)
        self.current_range = tk.Text(self.frame2, height=1, width=7)
        self.current_range.insert('1.0', '1e-6')
        self.current_range.grid(row=2, column=1)
        
        
        
        # Create waveform from result
        self.make_waveform_button = tk.Button(self.frame2, 
                                              text='Create waveform from last measurement', 
                                              command=self.make_waveform)
        self.make_waveform_button.grid(row=3, column=0, columnspan=2)
        
        
        
        # File directory
        path = os.path.abspath(os.path.dirname(sys.argv[0]))
        def_path = os.path.join(path, '\\rigol_waves\\')
        text = tk.Label(self.frame2, text='CSV file directory:')
        text.grid(row=4, column = 0)
        self.csv_dir = tk.Text(self.frame2, height=1, width=7)
        self.csv_dir.insert('1.0', def_path)
        self.csv_dir.grid(row=5, column=0, columnspan=2)
        
        
        # Recording/ processing parameters
        '''        
        Show previous recorded signal button
        
        
        To add:
            
        Create new waveform from result
         
        Save button
        
        Scope determines vertical scaling to use dynamically
        
        '''
        
    def get_units(self, n):    
        if n >= 1e-6 and n < 1e-3:
            return ('u', 1e-6)
        
        if n >= 1e-3 and n < 0:
            return ('m', 1e-3)
        
        if n >= 0 and n <= 1000:
            return ('', 1)
        
        if n > 1e3 and n <= 1e6:
            return ('k', 1e3)
        
        if n > 1e6:
            return ('M', 1e6)    
    
    
    
    def make_waveform(self):
        if not self.ft:
            print('No previous scan. Record data then try again')
            return
        
        # Get current waveform's phases from f_.csv file
        waveform = self.ft[0].waveform
        phase_file = 'f' + waveform[10:]
        phase_file = os.path.join(rigol_waves, phase_file)
        
        df = pd.read_csv(phase_file, skiprows=1, names=('index', 'f', 'phase'))
        phases = df['phase'].to_numpy()
        
        # Get |Z| and f from previous data set
        Z     = self.ft[0].Z
        freqs = self.ft[0].freqs
        
        # Average Z over all frames
        for i in self.ft:
            Z = np.mean(np.array([Z, self.ft[i].Z]), axis=0)
            
        amps = np.absolute(Z)
        
        # Create new waveform
        S, fname = create_waveform.Rigol_waveform(freqs, phases, 
                                    sample_freq=100000, total_time=1, 
                                    amax = 1, amps= amps, save_path = rigol_waves)
        
        # Display the new waveform
        self.show_waveform(new_waveform = S)
        
        
        # Reinitialize file selection list to incude new file
        self.file_list = [file for file in os.listdir(rigol_waves) 
                          if file.endswith('freqs.csv')
                          if file.startswith('Rigol')]
        
        del self.waveform_selector
        self.waveform_selector = tk.OptionMenu(self.frame, self.waveform, 
                                               *self.file_list, command=self.show_waveform)
        self.waveform_selector.grid(row=2, column=2)
    
    
    
    def show_waveform(self, selection=None, new_waveform=None):
        # Plot currently selected waveform
        try:
            if len(new_waveform) > 0:
                s = new_waveform
        
        except:
            file = os.path.join(rigol_waves, self.waveform.get())
            df = pd.read_csv(file, skiprows=9, names=('x', 'V'))
            s = df['V'].to_numpy()
        
        self.ax.clear()
        self.ax.plot(100000*np.fft.rfftfreq(len(s)), 
                     np.abs(np.fft.rfft(s)))
        self.ax.set_xlabel('Frequency/ Hz')
        self.ax.set_ylabel('Amplitude')
        self.ax.set_xscale('log')
        self.fig.tight_layout()
        self.canvas.draw_idle()
        
        
    
    
    def apply_waveform(self):
        # Send currently selected waveform to Arb
        file = os.path.join(rigol_waves, self.waveform.get())    
        df = pd.read_csv(file, skiprows=9, names=('x', 'V'))    
        s = df['V'].to_numpy()
        
        inst = self.rm.open_resource(self.arb.get())
        
        Vpp = self.waveform_vpp.get('1.0', 'end')
        Vpp = str(float(Vpp)/1000)
        
        try:
            inst.write('*RST')
            rigol_control.apply_waveform(inst, 1, s, Vpp)
        except:
            pass

        

        
        
    

    def record_signals(self):
        # Get recording time
        t = self.recording_time.get('1.0', 'end')
        current_range = self.current_range.get('1.0', 'end')
        current_range = float(current_range)
        
        
        
        try:
            t = float(t)
            t > 0
        except:
            print('Invalid time. Must be a real number > 0.')
            return
        
        # Initialize graph
        self.ax.clear()
        self.ax.set_xscale('log')
        self.ax.set_xlabel('Frequency/ Hz')
        self.ax.set_ylabel('Z/ $\Omega$')
        line1, = self.ax.plot([],[], 'o-')
        self.fig.tight_layout()
        self.canvas.draw_idle()
        
        # Connect to scope
        inst = self.rm.open_resource(self.scope.get())
             
        
        # Set and record some scope parameters
        inst.write('TRMD AUTO')
        inst.write('MSIZ 70K')
        inst.write('TDIV 100MS')
        inst.write('TRMD STOP')
        
        vdiv1       = float(inst.query('C1:VDIV?')[8:-2])
        voffset1    = float(inst.query('C1:OFST?')[8:-2])
        vdiv2       = float(inst.query('C2:VDIV?')[8:-2])
        voffset2    = float(inst.query('C2:OFST?')[8:-2])
        sara        = float(inst.query('SARA?')[5:-5])
        tdiv        = float(inst.query('TDIV?')[5:-2])
        frame_time  = 14*tdiv
        
        
        # Get applied frequencies
        file  = os.path.join(rigol_waves, self.waveform.get())    
        df    = pd.read_csv(file, skiprows=9, names=('x', 'V'))    
        
        freqs = 100000*np.fft.rfftfreq(len(df))
        V     = np.fft.rfft(df['V'])
        ftdf  = pd.DataFrame({
                'freqs': freqs,
                'V': V})
        
        ftdf  = ftdf[np.abs(ftdf['V']) > 100]
        applied_freqs = ftdf['freqs'].to_numpy()
            
            
        # Record starting time
        start_time = time.time()
        self.ft = {}
        frame = 0
        
        
        # Record frames
        print('')
        print('Recording for ~%d s' %t)
        while time.time() - start_time < t:
            d = siglent_control.record_single(inst, start_time, frame_time,
                                          vdiv1, voffset1, vdiv2, voffset2,
                                          sara, sample_time=1)
            print(f'Frame %s: {d.time:.2f} s'%frame)
            
            V = d.CH1data
            I = -d.CH2data * current_range
            
            Z = V/I
            
            df = pd.DataFrame({
                'freqs': d.freqs,
                'Z': Z})
            
            df = df[df['freqs'].isin(applied_freqs)]
            d.freqs = df['freqs'].to_numpy()
            d.Z = df['Z'].to_numpy()
            
            ydata = np.abs(d.Z)
            
            line1.set_xdata(d.freqs)          
            line1.set_ydata(ydata)
            
            self.ax.set_xlim(0.7*min(d.freqs), 1.3*max(d.freqs))
            self.ax.set_ylim(min(ydata)-1.05*min(ydata), 1.05*max(ydata))
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
            
            d.waveform = self.waveform.get()
            self.ft[frame] = d
            frame += 1
        
        print(f'Measurement complete. Total time {time.time()-start_time:.2f} s')
        

        
        


root = tk.Tk()
gui = MainWindow(root)
root.mainloop()
