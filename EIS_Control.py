import tkinter as tk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import sys
import time
import cmath
from datetime import date, datetime
from array import array
import pyvisa
from EIS_Control import rigol_control, siglent_control, create_waveform
from EIS_Control.funcs.recording_inits import init_recording, init_save
from EIS_Control.funcs.recording_funcs import FourierTransformData, record_frame, process_frame, save_frame
from EIS_Fit import EIS_fit
# from siglent_control import FourierTransformData
default_stdout = sys.stdout
default_stdin  = sys.stdin
default_stderr = sys.stderr



this_dir = rigol_control.__file__[:-16]
config_file = os.path.join(this_dir, 'config')


# Find dir with waveforms
rigol_waves = os.path.join(this_dir, 'waveforms')


# Get matplotlib style sheet and color cycle
plt.style.use(os.path.join(this_dir[:-13], 'scientific.mplstyle'))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


'''
To add:


'''



##### Create log file #####

LOGGING = True

log_file = os.path.expanduser('~\Desktop\EIS Output\log.txt')

#def log(text):
#    global log_file
#    if text != '\n' and text != '' and LOGGING:
#        with open(log_file, 'a') as f:
#            t = str(datetime.now().time())
#            f.write(t + '\t' + text + '\n')
#            f.close()

            
            
##### PrintLogger class #####

class PrintLogger(): 
    # Class to print console output into Tkinter window
    def __init__(self, Rec, textbox): # pass reference to text widget
        self.Rec = Rec
        self.textbox = textbox # keep ref

    def write(self, text):
        # use Recorder.log() to log console output
        self.Rec.log(text)
        self.textbox.insert(tk.END, text) # write text to textbox
        self.textbox.see('end') # scroll to end

    def flush(self): # needed for file like object
        pass



def TKObject(obj:str, frame, pos:tuple,
             variable=None, default=None, options=None,
             command=None, text=None, columnspan=1,
             size=None, returnselector=False):
    # Wrapper function for creating TKinter objects
    if obj == "Label":
        tkobj = tk.Label(frame, text=text)
        tkobj.grid(row=pos[0], column=pos[1], columnspan=columnspan)
    
    elif obj == 'StringVar':
        if not default:
            default = options[0]
        variable = tk.StringVar(frame)
        variable.set(default)
        selector = tk.OptionMenu(frame, variable, *options, command=command)
        selector.grid(row=pos[0], column=pos[1], columnspan=columnspan)
        if returnselector:
            return variable, selector
        else:
            return variable
        
    elif obj == 'Button':
        tkobj = tk.Button(frame, text=text, command=command)
        tkobj.grid(row=pos[0], column=pos[1], columnspan=columnspan)
        
    elif obj == 'IntVar':
        variable = tk.IntVar(value=default)
        option = tk.Checkbutton(frame, text=text,
                                variable=variable)
        option.grid(row=pos[0], column=pos[1], columnspan=columnspan)
        return variable
    
    elif obj == 'Text':
       variable = tk.Text(frame, height=size[0], width=size[1])
       variable.insert('1.0', default)
       variable.grid(row=pos[0], column=pos[1])
       return variable
   
    else:
        print(f'Unsupported obj: {obj}')




##### Recorder class #####

class Recorder:
    
    global this_dir, rigol_waves
    
    def __init__(self, root):
        self.root = root
        root.title("Arb/ Scope control")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    

        self.rm = pyvisa.ResourceManager()
        self.ft = None
        self.config = {}
        self.load_config()
        
        
        self.last_file_name = ''
        self.test_mode = False
        
        
        self.potentiostat   = None
        self.waveform       = None
        self.arb            = None
        self.scope          = None
        self.plot_Z         = None
        self.plot_phase     = None
        self.waveform_vpp   = None
        self.recording_time = None
        self.current_range  = None
        self.DC_offset      = None
        self.ref_corr_val   = None
        self.ref_corr_var   = None
        self.fit            = None
        self.circuit        = None
        self.time_plot_param= None
        self.time_plot_val  = None
            
        
        # Initialize frames and canvas
        
        # frame: upper left
        self.frame = tk.Frame(self.root)
        self.frame.grid(row=0, column = 0)
        
        # frame2: upper middle
        self.frame2 = tk.Frame(self.root)
        self.frame2.grid(row=0, column = 1)
        
        # frame3: lower middle (console output)
        self.frame3 = tk.Frame(self.root)
        self.frame3.grid(row=1, column=1)
        self.frame3.pack_propagate(0)
        
        # frame4: upper right
        self.frame4 = tk.Frame(self.root)
        self.frame4.grid(row=0, column=2)
        self.frame4.pack_propagate(0)
        
        # console printout to frame3
        self.console = tk.Text(self.frame3, width=50, height=25)
        self.console.grid(row=0, column=0)
        pl = PrintLogger(self, self.console)
        sys.stdout = pl
       
        # fig: lower left
        self.fig = plt.Figure(figsize=(5,4), dpi=100)
        self.ax  = self.fig.add_subplot(111)
        self.ax2 = self.ax.twinx()
        self.ax2.set_yticks([])
                
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=1, column=0)
        
        self.init_time_plot(1)
                
        
        
        #################################################
        ### FRAME 1: Instrument selection and control ###
        #################################################
        
        # Potentiostat selection
        
        TKObject('Label', self.frame,
                 text='Potentiostat:', pos=(1,1))
        self.potentiostat = TKObject('StringVar', self.frame,
                                      variable=self.potentiostat,
                                      default=self.config['potentiostat'],
                                      options=['Gamry', 'Autolab'], pos=(1,2))
                
        
        # Update config file button
        TKObject('Button', self.frame, text='Update config',
                 command=self.update_config_file, pos=(1,3))
        
        
        
        # Waveform selection dropdown menu
        self.file_list = [file for file in os.listdir(rigol_waves) 
                          if file.endswith('freqs.csv')
                          if file.startswith('Rigol')]
        TKObject('Label', self.frame,
                  text='Waveform:', pos=(2,1))
        
        self.waveform = TKObject('StringVar', self.frame,
                                  default=self.config['waveform'],
                                  variable=self.waveform, options=self.file_list,
                                  command=self.show_waveform, pos=(2,2))
        

                
        
        # VISA selection menu: Arb
        try:
            # Look for Rigol arb, i.e.
            # 'USB0::0x1AB1::0x0643::DG8A232302748::INSTR'
            default_arb = [inst for inst in self.rm.list_resources() 
                            if len(inst.split('::')) > 3 
                            and inst.split('::')[3].startswith('DG')][0]
        
        except:
            default_arb = ''
        
        try:
            TKObject('Label', self.frame, text='Arb:', pos=(3,1))
            self.arb = TKObject('StringVar', self.frame, default=default_arb,
                                      variable=self.arb, options=self.rm.list_resources(),
                                      pos=(3,2))
            TKObject('Button', self.frame, text='Apply Wave', command=self.apply_waveform, pos=(2,3))
        except:
             # If no instrument connected, no button
            self.test_mode = True
            pass
                
        
        
        # VISA selection menus: Scope
        try:
            # Look for Siglent scope
            default_scope = [inst for inst in self.rm.list_resources() 
                                if len(inst.split('::')) > 3 
                                and inst.split('::')[3].startswith('SDS')][0]
        
        except:
            default_scope = ''
        
        try:
            
            TKObject('Label', self.frame, text='Scope:', pos=(4,1))
            self.scope = TKObject('StringVar', self.frame, default=default_scope,
                                      variable=self.scope, options=self.rm.list_resources(),
                                      pos=(4,2))
        except:
             # If no instrument connected, no button
            self.test_mode = True
            pass
       
        
        
        
        # Record, save buttons
        TKObject('Button', self.frame, text='Record Signals', 
                 command=self.record_signals, pos=(3,3))
        TKObject('Button', self.frame, text='Record Reference', 
                 command=self.record_reference, pos=(5,1))
        TKObject('Button', self.frame, text='Save last measurement', 
                 command=self.save_last, pos=(5,2))
        TKObject('Button', self.frame, text='Multiplex', 
                 command=self.multiplex, pos=(5,3))
        TKObject('Button', self.frame, text='Record and save', 
                 command=self.record_and_save, pos=(4,3))
                
        
        self.plot_Z = TKObject('IntVar', self.frame, text='|Z|',
                               default=self.config['plot_Z'],
                               variable=self.plot_Z, pos=(6,1))
        
        
        self.plot_phase = TKObject('IntVar', self.frame, text='Phase',
                               default=self.config['plot_phase'],
                               variable=self.plot_phase, pos=(6,2))
                
        
        
        ########################################
        ### FRAME 2: User-adjustable options ###
        ########################################
        

        # Applied waveform amplitude
        TKObject('Label', self.frame2, text='Waveform Vpp (mV):', pos=(0,0))
        self.waveform_vpp = TKObject('Text', self.frame2, size=(1,7), pos=(0,1),
                                     default=self.config['waveform_vpp'])

        
        # Recording duration
        TKObject('Label', self.frame2, text='Recording time (s):', pos=(1,0))
        self.recording_time = TKObject('Text', self.frame2, size=(1,7), pos=(1,1),
                                     default=self.config['recording_time'])

        
        # Potentiostat current range
        TKObject('Label', self.frame2, text='Current range:', pos=(2,0))
        self.current_range = TKObject('Text', self.frame2, size=(1,7), pos=(2,1),
                                     default=self.config['current_range'])

        
        # DC Voltage offset
        TKObject('Label', self.frame2, text='DC Voltage (V):', pos=(3,0))
        self.DC_offset = TKObject('Text', self.frame2, size=(1,7), pos=(3,1),
                                     default=self.config['DC_offset'])
        
        
        TKObject('Button', self.frame2, text='Apply offset', pos=(3,2),
                 command = self.apply_offset)
         
        
        # Apply calibration correction
        TKObject('Label', self.frame2, text='Apply reference correction:', pos=(4,0))
        self.ref_corr_val = TKObject('Text', self.frame2, size=(1,7), pos=(4,1),
                                     default=self.config['ref_corr_val'])
        self.ref_corr_var = TKObject('IntVar', self.frame2, pos=(4,2),
                                     variable=self.ref_corr_var,
                                     default=self.config['ref_corr_var'])
                
        
        # Create waveform from result
        TKObject('Button', self.frame2, pos=(5,0), columnspan=2,
                 text='Create waveform from last measurement', 
                 command=self.make_waveform)       
        
        
        # Save options
        self.fit = TKObject('IntVar', self.frame2, pos=(6,0), text='Fit',
                                     variable=self.fit,
                                     default=self.config['fit'])
                
        
        self.circuit = TKObject('StringVar', self.frame2, variable=self.circuit,
                                options=['RRC', 'Randles_adsorption'],
                                pos=(6,1), columnspan=2)
                
        
        ########################################
        ###      X VS TIME PLOTTING OPTIONS  ###
        ########################################
        
        
        self.time_plot_param = TKObject('StringVar', self.frame4, pos=(0,2),
                                        variable=self.time_plot_param,
                                        default=self.config['time_plot_param'],
                                        options=['Phase', '|Z|', 'Parameter'],
                                        command=self.selector_changed)
        
        self.time_plot_val, self.time_plot_val_selector = TKObject('StringVar', 
                                      self.frame4, pos=(0,3),
                                      variable=self.time_plot_val,
                                      options=['-'], returnselector=True)
        
        
        self.selector_changed()
        self.update_config_file()
        
        ### END __INIT__ ###
    
    def log(self, text):
        global log_file
        if text != '\n' and text != '' and LOGGING:
            with open(log_file, 'a') as f:
                t = str(datetime.now().time())
                f.write(t + '\t' + text + '\n')
                f.close()    
    
        
    def load_config(self):
        if not os.path.isfile(config_file):
            print('No config file, using default settings')
            # Default options
            self.config = {
             'last_file_name': '',
             'test_mode': False,
             'potentiostat': 'Autolab',
             'waveform': 'Rigol_100k_1k_1_16freqs.csv',
             'plot_Z': 1,
             'plot_phase': 1,
             'waveform_vpp': '20',
             'recording_time': '10',
             'current_range': '1e-6',
             'DC_offset': '0.0',
             'ref_corr_val': '10k',
             'ref_corr_var': 1,
             'fit': 0,
             'circuit': 'RRC',
             'time_plot_param': 'Phase',
             'time_plot_val': '1'
             }
            
        else:
            print("Loaded config file")
            with open(config_file, 'r') as f:
                for line in f:              
                    key, val = line.split(':')
                    
                    if key == 'recording_time' and float(val) < 10:
                        val = '10'
                    
                    self.config[key] = val.strip('\n')
            
    
    
                
    def update_config_file(self):
        self.config = {'last_file_name': self.last_file_name,
                 'test_mode': self.test_mode,
                 'potentiostat': self.potentiostat.get(),
                 'waveform': self.waveform.get(),
                 'plot_Z': self.plot_Z.get(),
                 'plot_phase': self.plot_phase.get(),
                 'waveform_vpp': self.waveform_vpp.get('1.0', 'end').strip('\n'),
                 'recording_time': self.recording_time.get('1.0', 'end').strip('\n'),
                 'current_range': self.current_range.get('1.0', 'end').strip('\n'),
                 'DC_offset': self.DC_offset.get('1.0', 'end').strip('\n'),
                 'ref_corr_val': self.ref_corr_val.get('1.0', 'end').strip('\n'),
                 'ref_corr_var': self.ref_corr_var.get(),
                 'fit': self.fit.get(),
                 'circuit': self.circuit.get(),
                 'time_plot_param': self.time_plot_param.get(),
                 'time_plot_val': self.time_plot_val.get()
                 }
        
        with open(config_file, 'w') as f:
            for key, val in self.config.items():
                f.write(str(key) + ':' + str(val))
                f.write('\n')
            
            f.close()
        return
    
    
    
    def selector_changed(self):
        # Reinitialize parameter to plot vs time        
        self.get_waveform()  # Initialize list of frequencies
        _, _, params = self.initialize_circuit() # Initialize list of circuit elements
        
        
        if self.time_plot_param.get() in ('Phase', '|Z|'):
            time_plot_val = self.freqs
            
        
        elif self.fit.get():
            time_plot_val = params
        else:
            time_plot_val = ['-']
        
        menu = self.time_plot_val_selector['menu']
        menu.delete(0, 'end')
        for string in time_plot_val:
            menu.add_command(label=string, 
                             command=lambda value=string: 
                                 self.time_plot_val.set(value))
                
        
        self.time_plot_val.set(time_plot_val[0])
        self.update_config_file()
        
        
        
    def get_waveform(self):
        file = os.path.join(rigol_waves, self.waveform.get())
        df = pd.read_csv(file, skiprows=9, names=('x', 'V'))
        s = df['V'].to_numpy()
        
        
        ftdf = pd.DataFrame({'freqs': 100000*np.fft.rfftfreq(len(s)),
                             'V': np.abs(np.fft.rfft(s))})
        
        ftdf  = ftdf[np.abs(ftdf['V']) > 100]
        applied_freqs = ftdf['freqs'].to_numpy()
        
        self.freqs = applied_freqs.astype(int)
        self.update_config_file()
        return [s, applied_freqs]
    
    
    def get_freqs(self):
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
        
        return applied_freqs
    
    
    def get_correction_values(self):
        # Path
        ref_dir = os.path.join(this_dir, 'reference waveforms\\')
        
        # Resistance
        R = self.ref_corr_val.get('1.0', 'end')
        R = R[:-1]
        
        # Waveform
        waveform = self.waveform.get()
        
        if waveform.split('_')[1] == 'opt':
            # Use same correction factors for optimized waveform as
            # for unoptimized
            waveform = waveform.split('_')
            del waveform[1]
            waveform = '_'.join(waveform)
        
        
        # Put it all together and get corrections
        fname = 'REF_%s_%s'%(R, waveform)
        file = os.path.join(ref_dir, fname)
        
        try:
            corr_df = pd.read_csv(file, skiprows=1, names=('freqs', 
                                                           'Z_corr', 
                                                           'phase_corr')
                                  )
            
            Z_corr = corr_df['Z_corr'].to_numpy()
            phase_corr = corr_df['phase_corr'].to_numpy()
            return Z_corr, phase_corr
        
        except:
            print('Invalid reference file: ')
            print(file)
            print('Uncheck "Apply reference correction" or record a')
            print('reference spectrum of a resistor.\n')
            return 0, 0
        
    
      
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
        
        # Get f from previous data set
        # Average Z over all frames
        Z     = np.mean([self.ft[i].Z for i in self.ft], axis=0)
        freqs = self.ft[0].freqs
            
        amps = np.sqrt(np.absolute(Z))
        
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
        self.update_config_file()
        print('\n')
        
    
    
    def show_waveform(self, selection=None, new_waveform=None):
        # Plot currently selected waveform
        try:
            if len(new_waveform) > 0:
                s = new_waveform
        
        except:
            [s, _] = self.get_waveform()
        
        self.ax.set_xscale('linear')
        self.ax.clear()
        self.ax2.clear()
        self.ax2.set_yticks([])
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
        signal = df['V'].to_numpy()
        
        inst = self.rm.open_resource(self.arb.get())
        
        Vpp = self.waveform_vpp.get('1.0', 'end')
        Vpp = str(float(Vpp)*2/1000)
        
        try:
            inst.write('*RST')
            rigol_control.apply_waveform(inst, 1, signal, Vpp)
        except:
            print('Could not connect')
            pass
        
        self.update_config_file()
    
    
    
    def apply_offset(self):
        
        try:
            # Connect to scope
            inst = self.rm.open_resource(self.scope.get())
                 
            # Set scope parameters
            # inst.write('C1:VDIV 5mV')
            inst.write('C1:OFST %s' %self.DC_offset.get('1.0', 'end'))
            
            inst.write('TRMD AUTO')
            
        
        except:
            # No instrument connected
            pass
        
        self.update_config_file()
    
    

    def init_time_plot(self, n_plots):
               
        if n_plots == 1:
            self.timefig = plt.Figure(figsize=(5,4), dpi=100)
            self.timeax = self.timefig.add_subplot(111)
            
            self.timeax.plot([],[], 'ok')
            self.timeax.set_xlabel('Time/ s')
            
        elif n_plots > 1:
            # Calling plt.subplots() directly makes figures print to
            # console when GUI is closed for some reason... this way doesn't
            self.timefig = plt.Figure(figsize=(5,4), dpi=100)
            self.timeax = []
            for i in range(n_plots):
                self.timeax.append(self.timefig.add_subplot(n_plots, 1, i+1))
            self.timeax = np.array(self.timeax)
            
            self.timeax[-1].set_xlabel('Time/s')
            for ax in self.timeax:
                ax.plot([],[],'ok')
                
        
        if hasattr(self, 'timecanvas'):
            self.timecanvas.get_tk_widget().destroy()
                        
        self.timecanvas = FigureCanvasTkAgg(self.timefig, master=root)
        self.timecanvas.get_tk_widget().grid(row=1, column=2)
        
        self.timefig.tight_layout()
        self.timecanvas.draw()
        self.timecanvas.flush_events()
    
    
    
    def update_time_plot(self, t, freqs, Z, phase, params, ax=None):
        # Plot the most recent point of whatever's selected
        # (self.time_plot_param) onto the ... vs t graph
        
        if not ax:
            ax = self.timeax
        
        # Get the value
        if self.time_plot_param.get() == 'Phase':
            freq = int(self.time_plot_val.get())
            idx = np.where(freqs == freq)[0][0]
            val = phase[idx]
        
        elif self.time_plot_param.get() == '|Z|':
            freq = int(self.time_plot_val.get())
            idx = np.where(freqs == freq)[0][0]
            val = Z[idx]
        
        elif self.time_plot_param.get() == '-':
            return
        
        else:
            freq = None
            param = self.time_plot_val.get()
            val = params[param]
        
        
        # Plot the value
        if freq:
            label = f'{self.time_plot_param.get()} @ {freq} Hz'
        else:
            label = f'{self.time_plot_val.get()}'
        
        ax.plot(t, val, 'ok')
        ax.set_ylabel(label)
#        self.timefig.tight_layout()
        self.timefig.canvas.draw_idle()
#        self.timefig.canvas.flush_events()
            
            
    
    def initialize_circuit(self):
        # Returns current EEC parameters, bounds, and starting guess for fit
        # Update as more circuits are added
        
        if self.circuit.get() == 'Randles_adsorption':
                bounds = {
                        'R1': [1e-1, 1e9],
                        'R2': [1e-1, 1e9],
                        'Q1': [1e-15, 1],
                        'n1': [0.9,1.1],
                        'Q2': [1e-15, 1],
                        'n2': [0.8,1.1]
                        }
                
                starting_guess = {
                    'R1': 284, 
                    'R2': 62000, 
                    'Q1': 1.8e-07, 
                    'n1': 1, 
                    'Q2': 3.2e-07, 
                    'n2': 0.9
                    }
            
                    
        elif self.circuit.get() == 'RRC':
            bounds = {
                    'R1': [1e-1, 1e9],
                    'R2': [1e-1, 1e9],
                    'C': [1e-15, 1]
                    }
            
            starting_guess = None
          
        
        params = [param for param, val in bounds.items()]
        
        self.update_config_file()
        
        return bounds, starting_guess, params
        
        
        
        


    def record_signals(self, save=False, silent=True, plot_time_plot=True,
                       axes = None, new_time_plot=True, n_plots=1, **kwargs):
        '''
        Record impedance for one, non-multiplexed electrode for the set
        amount of time.

        Parameters
        ----------
        save : Bool, optional
            Whether or not to save data as ASCII. The default is False.
        silent : Bool, optional
            If True, don't printout to console. The default is True.

        Returns
        -------
        None.

        '''
        
        self.update_config_file()
        
        plot_Z     = self.plot_Z.get()
        plot_phase = self.plot_phase.get()
        
        
        # Get waveform correction factors
        if self.ref_corr_var.get():
            try:
                Z_corr, phase_corr = self.get_correction_values()
            except:
                # get_correction_values() returns 0 if 
                # invalid reference file
                return
                
        # Connect to scope
        inst = self.rm.open_resource(self.scope.get())
        
        # Initialize recording
        recording_params = init_recording(self, new_time_plot, 
                                          n_plots, save)
            
        # Get recording time
        t = self.recording_time.get('1.0', 'end')             
        
        try:
            t = float(t)
            t > 0
        except:
            print('Invalid time. Must be a real number > 0.')
            return
        
        
        
        # Initialize save files
        if save:
            # Make save folder
            try:
                name = tk.simpledialog.askstring('Save name', 'Input save name:',
                                             initialvalue = self.last_file_name)
                today = str(date.today())
                save_path = os.path.join(os.path.expanduser('~\Desktop\EIS Output'), 
                                       today, name)
                createFolder(save_path)
                self.last_file_name = name
            
            except:
                # User hits cancel
                return
            

            recording_files = init_save(self, save_path)
        
        else:
            recording_files = {}
            
#        recording_files = {'time_file':time_file,
#                           'meta_file':meta_file,
#                           'fits_file':fits_file,
#                           'DC_file':DC_file,
#                           'save_path':save_path}          
        
        ### RECORDING MAIN LOOP ###
        
        # Record starting time
        start_time = time.time()
        self.ft = {}
        
        # Record frames
        if not silent:
            print('')
            print('Recording for ~%d s' %t)
            
        frame = 0
        while time.time() - start_time < t:
            self.ft[frame] = record_frame(self, inst, 1.4*1.2, recording_params,
                                          time.time() - start_time, recording_files,
                                          frame, save=save)   
            if not silent:
                print(f'Frame {frame}: {self.ft[frame].time:.2f} s')                
            frame += 1
        
        inst.write('TRMD AUTO')
        
        # Process the final frame
        t, freqs, Z, phase, fits = process_frame(self, frame-1, 
                                                 update_time_plot=True)
           
        # if plot_time_plot:
        #     self.update_time_plot(t, freqs, Z, phase, fits)
        
        try:
            if not silent:
                print(self.ft[frame-1].params)
        except:
            pass
        
        print(f'Measurement complete. Total time {time.time()-start_time:.2f} s')
        
        if save:
            # Save last frame
            save_frame(self, frame-1, self.ft[frame-1], recording_files)
            
            # Save metadata
            self.fig.savefig(save_path+'\\0000_fig', dpi=100)
            
            with open(recording_files['meta_file'], 'a') as f:
                avg_Vpp = np.mean([self.ft[frame].Vpp for frame in self.ft])
                f.write(f'\nExperimental Vpp (V): {avg_Vpp}')
            
            print('Saved as ASCII:', save_path, '\n')
        
        
        
            
    def multiplex(self):
                        
        # Ask for titration/ invivo experiment
        exp_type = tk.simpledialog.askstring(
                        'Experiment type', 
                        'Titration (0) or in-vivo (1)?') 
        if exp_type == '0':
            exp_type = 'titration'
#            if int(self.recording_time.get('1.0', 'end')) != 10:
#                print('Set recording time to 10 for multiplexing!')
#                return
        elif exp_type == '1':
            exp_type = 'invivo'
        else:
            print('Invalid entry')
            return
        
        
        # Ask for number of multiplexed channels
        no_of_channels = tk.simpledialog.askstring(
                        'Number of channels', 'Number of channels:')    
        no_of_channels = int(no_of_channels)
        
        # Ask for electrode IDs (if desired, otherwise default to 1,2,3,..)
        elec_numbers = tk.simpledialog.askstring(
                        'Electrode numbers', 'Electrode numbers:',
                        initialvalue = ','.join([str(i) for i in range(1,no_of_channels+1)]))
        elec_numbers = elec_numbers.split(',')
        
        if exp_type == 'titration':
            # Ask for number of concentrations
            number_of_concs = tk.simpledialog.askstring(
                            'Number of concentrations', 'Number of concentrations:')
            number_of_concs = int(number_of_concs)
        
        if exp_type == 'invivo':
            recording_time = self.recording_time.get('1.0', 'end')
            recording_time = float(recording_time)
        
        # Path of triggering file
        updatefile = os.path.join(this_dir, 'update.txt')
        
        # Clean up trigger file to start
        if os.path.exists(updatefile):
            os.remove(updatefile)
        
                
        # Iterate through concentrations
        if exp_type == 'titration':
            for _ in range(number_of_concs):
                
                # Ask for concentration
                conc = tk.simpledialog.askstring(title=None,
                                                 prompt='Concentration: ')
                            
                # Wait for autolab to create start file
                # ULTRA bad way of triggering recording...
                while os.path.exists(updatefile) == False:
                    time.sleep(0.1)
                
                # Multiplex record and save
                for i in range(no_of_channels):
                    
                    while os.path.exists(updatefile) == False:
                        # Wait for autolab to create start file
                        self.root.after(5) #wait 5 ms
                        
                    
                    print(f'Recording electrode {elec_numbers[i]}, {conc}')
                    self.record_signals(silent=True)
                    self.save_last(name = f'{elec_numbers[i]}_{conc}')
                    
                    os.remove(updatefile)
        
        
        elif exp_type == 'invivo': 
            
            name = tk.simpledialog.askstring('Save name', 'Input save name:',
                                                 initialvalue = self.last_file_name)
            self.last_file_name = name
            
            conc = ''
            
            inst = self.rm.open_resource(self.scope.get())
            
            # self.recording_time.delete('1.0', 'end')
            # self.recording_time.insert('1.0', '1.5')
            
            # self.init_time_plot(no_of_channels)
            
            recording_params = init_recording(self, new_time_plot=True,
                                              n_plots=no_of_channels, save=True)
            
            
            
            while os.path.exists(updatefile) == False:
                self.root.after(5)
            
            # Start recording when Autolab makes update file
            s_t = time.strftime("%H%M%S", time.gmtime(time.time()))
            
            today = str(date.today())
            save_path = os.path.join(os.path.expanduser('~\Desktop\EIS Output'), 
                                       today, name)
            createFolder(save_path)
            
            recording_files = init_save(self, save_path)
            frame = 0            
              
            
            start_time = time.time()
            self.ft = {}
            self.log('Experiment starting')
            while time.time() - start_time < recording_time:
                # Multiplex
                for i in range(no_of_channels):
                    
                    
                    while os.path.exists(updatefile) == False:
                        # Wait for autolab to create start file
                        self.root.after(1)
                    
                    self.log(f'Starting electrode {elec_numbers[i]}, frame {frame}')
                    ftime = time.time() - start_time
                    ft = record_frame(self, inst, 1.4*1.2, recording_params,
                                      ftime, recording_files, frame, 
                                      save=True, process_last=True,
                                      ax=self.timeax[i], 
                                      multiplex_fname=f'{elec_numbers[i]}_{frame}')
                                    
                    
                    self.ft[frame] = ft
                    
                    frame += 1
#                    i += 1
                    
                    
                    
                    
                    # Record and save the frame
                    
                    # self.record_signals(silent=True, new_time_plot=False,
                    #                     savefig=False, plot_time_plot=False)
                    # self.save_last(name = f'{elec_numbers[i]}_{int(ftime)}',
                    #                subpath = s_t, savefig=False)
                    
                    # # Plot frame to time plot
                    # freqs = self.ft[frame].freqs
                    # Z     = self.ft[0].Z
                    # phase = self.ft[0].phase
                    # params= self.ft[0].params
                    # self.update_time_plot(ftime, freqs, Z, phase, params,
                    #                       ax = self.timeax[i])
                    
                    os.remove(updatefile)
            
            # Process and save the last frame
            process_frame(self, frame - 1, update_time_plot=True, 
                          ax=self.ft[frame-1].ax)
            save_frame(self, frame-1, self.ft[frame-1], recording_files,
                       multiplex_fname= self.ft[frame-1].name)
            print('Experiment finished')
#            self.recording_time.insert('1.0', str(recording_time))
                
            
                            
        
                
        
    def record_reference(self):
        # Record impedance spectrum of a resistor to calibrate
        # Save with resistance and waveform labelled
        
        # Check that reference correction is unchecked
        if self.ref_corr_var.get():
            print('Uncheck "Apply Reference Correction" before\n'+
                  'recording a reference spectrum!\n')  
                  
            return
        
        
        # Prompt for resistance value
        R = tk.simpledialog.askstring('Calibration', 'Resistance:')
        
        # Break on "cancel" button
        if not R:
            return
        
        
        # Record spectra
        self.record_signals()
        

        # Determine reference file path/ name
        ref_dir = os.path.join(this_dir, 'reference waveforms\\')
        
        waveform = self.waveform.get()
        
        name = 'REF_%s_%s'%(R, waveform)
        
        out_file = os.path.join(ref_dir, name)
        
        # Average spectra
        freqs = self.ft[1].freqs
        Z = np.mean([np.abs(self.ft[i].Z) for i in self.ft], axis=0)
        phase = np.mean([self.ft[i].phase for i in self.ft], axis=0)
        
        
        if R.endswith('k'):
            R = 1e3*float(R[:-1])
        
        elif R.endswith('M'):
            R = 1e6*float(R[:-1])
        
        else:
            R = float(R)
        
        
        # Determine corrections
        Z_corr      = Z / R
        phase_corr  = phase
                
        df = pd.DataFrame(
            {'freq': freqs,
            'Z_corr': Z_corr,
            'phase_corr': phase_corr}
            )
        
        # Save to csv
        df.to_csv(out_file, index=False)
        
        print('Saved correction file to:')
        print(out_file, '\n')
        
        
        
    def save_frame(self, num, freqs, re, im, save_path):
        d = pd.DataFrame(
            {'f': freqs,
            're': re,
            'im': im}
            )
        
        if num is None:
            num = 0
        
        try:
            fname = save_path + f'\\{num:04}s.txt'
        except:
            # Passing string as save name instead
            fname = save_path + f'\\{num}.txt'
    
        d.to_csv(fname, columns = ['f', 're', 'im'],
                     header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], 
                     sep = '\t', index = False, encoding='ascii')
        
        
        
    def save_last(self, name = None, savefig=True, subpath=''):
               
        if self.ft:
            try:
                if not name:
                    name = tk.simpledialog.askstring('Save name', 'Input save name:',
                                                 initialvalue = self.last_file_name)
                
                self.last_file_name = name
                
                today = str(date.today())
                                    
                folder_path = os.path.join(os.path.expanduser('~\Desktop\EIS Output'), 
                                           today, subpath, name)
                
                createFolder(folder_path)
                
                meta_file = os.path.join(folder_path, '0000_Metadata.txt')
                time_file = os.path.join(folder_path, '0000_time_list.txt')
                DC_file = os.path.join(folder_path, '0000_DC_currents.txt')
                
                
                with open(time_file, 'w') as f:
                    for i, _ in self.ft.items():
                        ftime = str(self.ft[i].time)
                        f.write(ftime + '\n')
                
                with open(meta_file, 'w') as f:
                    f.write('Waveform Vpp (mV): '+ str(self.waveform_vpp.get('1.0', 'end')))
                    f.write('Waveform: '+ str(self.waveform.get()))
                    s_t = datetime.now().strftime("%H:%M:%S.%f")
                    f.write('\nStart time: %s'%s_t)
                    avg_Vpp = np.mean([self.ft[frame].Vpp for frame in self.ft])
                    f.write(f'\nExperimental Vpp (V): {avg_Vpp}')
                    
                with open(DC_file, 'w') as f:
                    for i, _ in self.ft.items():
                        DC = self.ft[i].mean_I
                        f.write(DC + '\n')    
                  
                    
                for i, _ in self.ft.items():
                    re = np.real(self.ft[i].Z)
                    im = np.imag(self.ft[i].Z)
                    freqs = self.ft[i].freqs
                    
                    self.save_frame(i, freqs, re, im, folder_path)
                
                
                if savefig:    
                    self.fig.savefig(folder_path+'\\0000_fig', dpi=100)
                
                print('Saved as ASCII:', folder_path, '\n')
                 
            
            except:
                # User hits cancel
                # Still option to save previous run
                pass
                

        else:
            print('No previous measurement to export\n')


    

    def record_and_save(self):
        self.record_signals(save=True)
     
        
         
        
        
    ########################################
    ###     END OF MAINWINDOW CLASS      ###
    ########################################


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)


root = tk.Tk()
gui = Recorder(root)
root.mainloop()

sys.stdout = default_stdout
sys.stdin = default_stdin
sys.stderr = default_stderr


