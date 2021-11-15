import tkinter as tk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import sys
import time
from datetime import date
import pyvisa
import rigol_control
import siglent_control
import create_waveform
default_stdout = sys.stdout
default_stdin  = sys.stdin
default_stderr = sys.stderr


this_dir = rigol_control.__file__[:-16]


# Find dir with waveforms
rigol_waves = os.path.join(this_dir, 'waveforms')


# Get matplotlib style sheet and color cycle
plt.style.use(os.path.join(this_dir[:-13], 'scientific.mplstyle'))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


'''
To add:
    
Vdc offset

Record reference spectrum

'''

class PrintLogger(): 
    # Class to print console output into Tkinter window
    def __init__(self, textbox): # pass reference to text widget
        self.textbox = textbox # keep ref

    def write(self, text):
        self.textbox.insert(tk.END, text) # write text to textbox
        self.textbox.see('end') # scroll to end

    def flush(self): # needed for file like object
        pass



class MainWindow:
    
    global this_dir, rigol_waves
    
    def __init__(self, root):
        self.root = root
        root.title("Arb/ Scope control")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    

        self.rm = pyvisa.ResourceManager()
        self.ft = None
            
        
        # Initialize frames and canvas
        
        # frame: upper left
        self.frame = tk.Frame(self.root)
        self.frame.grid(row=0, column = 0)
        
        # frame2: upper right
        self.frame2 = tk.Frame(self.root)
        self.frame2.grid(row=0, column = 1)
        
        # frame3: lower right (console output)
        self.frame3 = tk.Frame(self.root)
        self.frame3.grid(row=1, column=1)
        self.frame3.pack_propagate(0)
        
        # console printout to frame3
        self.console = tk.Text(self.frame3, width=50, height=25)
        self.console.grid(row=0, column=0)
        pl = PrintLogger(self.console)
        sys.stdout = pl
       
        # fig: lower left
        self.fig = plt.Figure(figsize=(5,4), dpi=100)
        self.ax  = self.fig.add_subplot(111)
        self.ax2 = self.ax.twinx()
        self.ax2.set_yticks([])
                
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=1, column=0)
        
        # Other vars to initialize
        self.last_file_name = ''
        
        
        
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
        try:
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
            
        except:
            # If no instrument connected, no button
            self.test_mode = True
            pass
        
        
        
        # VISA selection menus: Scope
        try:
            text = tk.Label(self.frame, text='Scope:')
            text.grid(row=4, column=1)
            self.scope = tk.StringVar(self.frame)
            self.scope.set('USB0::0xF4ED::0xEE3A::SDS1EDEX5R5381::INSTR')
            self.scope_selector = tk.OptionMenu(self.frame, self.scope, 
                                                   *self.rm.list_resources())
            self.scope_selector.grid(row=4, column=2)
        except:
            pass
        
        
        
        # Record, save buttons
        self.record_signals_button = tk.Button(self.frame, text='Record Signals', 
                                               command=self.record_signals)
        self.record_signals_button.grid(row=4, column=3)
        
        
        self.save_button = tk.Button(self.frame, text='Save last measurement', 
                                                   command=self.save_last)
        self.save_button.grid(row=5, column=3, columnspan=2)
        

        self.record_save_button = tk.Button(self.frame, text='Record and save', 
                                                   command=self.record_and_save)
        self.record_save_button.grid(row=5, column=2, columnspan=2)
        
        
        
        # Plot Z, phase toggles
        self.plot_Z = tk.IntVar(value=1)
        self.plot_Z_option = tk.Checkbutton(self.frame, text='|Z|', 
                                                variable=self.plot_Z)
        self.plot_Z_option.grid(row=6, column=1)
        
        
        self.plot_phase = tk.IntVar(value=1)
        self.plot_phase_option = tk.Checkbutton(self.frame, text='Phase', 
                                                variable=self.plot_phase)
        self.plot_phase_option.grid(row=6, column=2)
        
        
        
        
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
        
        
        
        # Potentiostat current range
        text = tk.Label(self.frame2, text='Current range:')
        text.grid(row=2, column = 0)
        self.current_range = tk.Text(self.frame2, height=1, width=7)
        self.current_range.insert('1.0', '1e-6')
        self.current_range.grid(row=2, column=1)
        
        
        
        # DC Voltage offset
        text = tk.Label(self.frame2, text='DC Voltage:')
        text.grid(row=3, column = 0)
        self.DC_offset = tk.Text(self.frame2, height=1, width=7)
        self.DC_offset.insert('1.0', '0.0')
        self.DC_offset.grid(row=3, column=1)
        
        self.DC_offset_button = tk.Button(self.frame2, 
                                              text='Apply offset', 
                                              command=self.apply_offset)
        self.DC_offset_button.grid(row=3, column=2)
        
        
        
        # Create waveform from result
        self.make_waveform_button = tk.Button(self.frame2, 
                                              text='Create waveform from last measurement', 
                                              command=self.make_waveform)
        self.make_waveform_button.grid(row=4, column=0, columnspan=2)
        
        
        
        
        # Save options
        
        text = tk.Label(self.frame2, text='Save as...')
        text.grid(row=6, column = 0)        
        
        self.asciiVar = tk.IntVar(value=1)
        self.save_ascii_option = tk.Checkbutton(self.frame2, text='ASCII', 
                                                variable=self.asciiVar)
        self.save_ascii_option.grid(row=6, column=1)
        
        self.csvVar = tk.IntVar(value=0)
        self.save_csv_option = tk.Checkbutton(self.frame2, text='CSV', 
                                              variable=self.csvVar)
        self.save_csv_option.grid(row=6, column=2)
        
                
        
        
        
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
        print('\n')
    
    
    
    def show_waveform(self, selection=None, new_waveform=None):
        # Plot currently selected waveform
        try:
            if len(new_waveform) > 0:
                s = new_waveform
        
        except:
            file = os.path.join(rigol_waves, self.waveform.get())
            df = pd.read_csv(file, skiprows=9, names=('x', 'V'))
            s = df['V'].to_numpy()
        
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
        s = df['V'].to_numpy()
        
        inst = self.rm.open_resource(self.arb.get())
        
        Vpp = self.waveform_vpp.get('1.0', 'end')
        Vpp = str(float(Vpp)/1000)
        
        try:
            inst.write('*RST')
            rigol_control.apply_waveform(inst, 1, s, Vpp)
        except:
            print('Could not connect')
            pass
    
    
    
    def apply_offset(self):
        
        try:
            # Connect to scope
            inst = self.rm.open_resource(self.scope.get())
                 
            
            # Set scope parameters
            inst.write('TRMD AUTO')
            inst.write('MSIZ 70K')
            inst.write('TDIV 100MS')
            inst.write('TRMD STOP')
            inst.write('C1:OFST %sV' %self.DC_offset)
        
        except:
            # No instrument connected
            pass
        


    def record_signals(self):
        
        plot_Z     = self.plot_Z.get()
        plot_phase = self.plot_phase.get()
        
        self.ax.set_xscale('linear')
        self.ax.clear()
        self.ax.set_xscale('linear')
        self.ax2.clear()
        
        line1, = self.ax.plot([],[], '-', color=colors[0])
        line2, = self.ax2.plot([],[], 'o', color=colors[1])
               
        
        self.ax.set_xscale('log')
        self.ax.set_xlabel('Frequency/ Hz')
        self.fig.tight_layout()
        self.canvas.draw_idle()
        
        
        try:
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
            
            
            # Connect to scope
            inst = self.rm.open_resource(self.scope.get())
                 
            
            # Set and record some scope parameters
            inst.write('TRMD AUTO')
            inst.write('MSIZ 70K')
            inst.write('TDIV 100MS')
            inst.write('TRMD STOP')
            inst.write('C1:OFST %sV' %self.DC_offset)
            
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
                phase = np.angle(V/I, deg=True)
                
                df = pd.DataFrame(
                        {
                        'freqs': d.freqs,
                        'Z': Z,
                        'phase': phase
                        }
                )
                
                df = df[df['freqs'].isin(applied_freqs)]
                d.freqs = df['freqs'].to_numpy()
                d.Z = df['Z'].to_numpy()
                d.phase = np.angle(d.Z, deg=True)
                
                # Plot this result to figure canvas
                Z = np.abs(d.Z)
                phase = d.phase
                
                # Determine which plot to make
                if plot_Z:
                    if not plot_phase:
                        line1.set_xdata(d.freqs)          
                        line1.set_ydata(Z)
                        self.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
                        self.ax.set_ylabel('|Z|/ $\Omega$')
                        self.ax2.set_yticks([])
                
                if plot_phase:
                    if not plot_Z:
                        line2.set_xdata(d.freqs)
                        line2.set_ydata(phase)
                        self.ax.set_ylim(min(phase)-10, max(phase)+10)
                        self.ax.set_ylabel('Phase/ $\degree$')
                        self.ax2.set_yticks([])
                    
                if plot_Z and plot_phase:
                    line1.set_xdata(d.freqs)
                    line2.set_xdata(d.freqs)
                    line1.set_ydata(Z)
                    line2.set_ydata(phase)
                    self.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
                    self.ax2.set_ylim(min(phase)-10, max(phase)+10)
                    self.ax.set_ylabel('|Z|/ $\Omega$')
                    self.ax2.set_ylabel('Phase/ $\degree$')
                    
                    
                # Draw the plot
                self.fig.tight_layout()
                self.ax.set_xticks([1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
                self.ax.set_xlim(0.7*min(d.freqs), 1.5*max(d.freqs))
                self.fig.canvas.draw()
                self.fig.canvas.flush_events()
                
                # Save the FT data
                d.waveform = self.waveform.get()
                self.ft[frame] = d
                frame += 1
            
            print(f'Measurement complete. Total time {time.time()-start_time:.2f} s\n')
        
        
        except:
            self.test_mode_record()
            
            
            
            
    def test_mode_record(self):

        plot_Z     = self.plot_Z.get()
        plot_phase = self.plot_phase.get()
        
        self.ax.set_xscale('linear')
        self.ax.clear()
        self.ax.set_xscale('linear')
        self.ax2.clear()
        
        line1, = self.ax.plot([],[], '-', color=colors[0])
        line2, = self.ax2.plot([],[], 'o', color=colors[1])
               
        
        self.ax.set_xscale('log')
        self.ax.set_xlabel('Frequency/ Hz')
        self.fig.tight_layout()
        self.canvas.draw_idle()
        
        if self.test_mode:
            # Show fake data if no instrument connected
            
            # Make up data
            d1 = siglent_control.FourierTransformData(
                time = 1.2, 
                freqs = np.logspace(1,3, num=20), 
                CH1data = [], 
                CH2data = [],
                Z = np.linspace(1,1000, num = 20) + 1j*np.linspace(1,1000, num=20),
                phase = np.linspace(1,1000, num=20),
                waveform = self.waveform.get()
                )
            
            d2 = siglent_control.FourierTransformData(
                time = 2.7, 
                freqs = np.logspace(1,3, num=20), 
                CH1data = [], 
                CH2data = [],
                Z = np.linspace(1,2000, num = 20) + 1j*np.linspace(1,2000, num=20),
                phase = np.linspace(1,2000, num=20),
                waveform = self.waveform.get()
                )
            
            # Data to plot
            Z = np.abs(d1.Z)
            phase = np.angle(Z, deg=True)
              
            if plot_Z:
                if not plot_phase:
                    line1.set_xdata(d1.freqs)          
                    line1.set_ydata(Z)
                    self.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
                    self.ax.set_ylabel('|Z|/ $\Omega$')
                    self.ax2.set_yticks([])
            
            if plot_phase:
                if not plot_Z:
                    line2.set_xdata(d1.freqs)
                    line2.set_ydata(phase)
                    # self.ax.set_ylim(min(phase)-1.05*min(phase), 1.05*max(phase))
                    self.ax.set_ylabel('Phase/ $\degree$')
                    self.ax2.set_yticks([])
                
            if plot_Z and plot_phase:
                line1.set_xdata(d1.freqs)
                line2.set_xdata(d1.freqs)
                line1.set_ydata(Z)
                line2.set_ydata(phase)
                self.ax.set_ylim(min(Z)-1.05*min(Z), 1.05*max(Z))
                # self.ax2.set_ylim(min(phase)-1.05*min(phase), 1.05*max(phase))
                self.ax.set_ylabel('|Z|/ $\Omega$')
                self.ax2.set_ylabel('Phase/ $\degree$')
                
                
            self.fig.tight_layout()
            self.ax.set_xticks([1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
            self.ax.set_xlim(0.7*min(d1.freqs), 1.5*max(d1.freqs))
            
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
            
            # Data to save
            self.ft = {
                0: d1,
                1: d2
                }
        
        
        
        
    def save_last(self):
        
        def createFolder(directory):
            try:
                if not os.path.exists(directory):
                    os.makedirs(directory)
            except OSError:
                print ('Error: Creating directory. ' +  directory)
        
        
        if self.ft:
            try:
                name = tk.simpledialog.askstring('Save name', 'Input save name:',
                                                 initialvalue = self.last_file_name)
                
                self.last_file_name = name
                
                today = str(date.today())
                
                if self.asciiVar.get():            
                    
                    folder_path = os.path.join(os.path.expanduser('~\Desktop\EIS Output'), 
                                               today, name)
                    
                    createFolder(folder_path)
                    
                    time_file = os.path.join(folder_path, '0000_time_list.txt')
                    
                    with open(time_file, 'w') as f:
                        for i, _ in self.ft.items():
                            time = str(self.ft[i].time)
                            f.write(time + '\n')
                        
                    f.close()
                        
                    for i, _ in self.ft.items():
                        re = np.real(self.ft[i].Z)
                        im = np.imag(self.ft[i].Z)
                        freqs = self.ft[i].freqs
                        
                        d = pd.DataFrame(
                            {'f': freqs,
                            're': re,
                            'im': im}
                            )
                        
                        fname = folder_path + '\\' + f'{i:04}' +'s.txt'
                    
                        d.to_csv(fname, columns = ['f', 're', 'im'],
                                     header = ['<Frequency>', '<Re(Z)>', '<Im(Z)>'], 
                                     sep = '\t', index = False, encoding='ascii')
                    
                    meta_file = os.path.join(folder_path, '0000_Metadata.txt')
                    
                    with open(meta_file, 'w') as f:
                        f.write('Waveform Vpp (mV): '+ str(self.waveform_vpp.get('1.0', 'end')))
                        f.write('Waveform: '+ str(self.waveform.get()))
                        
                    f.close()
                    
                        
                    self.fig.savefig(folder_path+'\\0000_fig', dpi=100)
                    
                    print('Saved as ASCII:', folder_path, '\n')
                 
                    
                if self.csvVar.get():
                    print('csv saving not yet supported...\n')
            
            except:
                # User hits cancel
                # Still option to save previous run
                pass
                

        else:
            print('No previous measurement to export\n')


    

    def record_and_save(self):
        self.record_signals()
        self.save_last()
            
        
        
        
    ########################################
    ###     END OF MAINWINDOW CLASS      ###
    ########################################
        

root = tk.Tk()
gui = MainWindow(root)
root.mainloop()
sys.stdout = default_stdout
sys.stdin = default_stdin
sys.stderr = default_stderr

