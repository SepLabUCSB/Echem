import tkinter as tk
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import pyvisa
import rigol_control
import siglent_control

csv_dir = r'C:\Users\BRoehrich\Desktop\git\echem\EIS_Control\csv'
rigol_waves = r'C:\Users\BRoehrich\Desktop\git\echem\EIS_Control\rigol_waves'



class MainWindow:
    
    global csv_dir, rigol_waves
    
    def __init__(self, root):
        self.root = root
        root.title("Arb/ Scope control")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    
        self.rm = pyvisa.ResourceManager()
        
        
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
        
        
        
        # Waveform selection dropdown menu
        text = tk.Label(self.frame, text='Waveform:')
        text.grid(row=2, column=1)
        self.file_list = [file for file in os.listdir(rigol_waves) 
                          if file.endswith('freqs.csv')]
        
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
        
        
        
        # Recording/ processing parameters
        '''
        Recording time (total_time)
        Waveform amplitude
        
        Current range for AD conversion
        
        
        To add:
            
        Create new waveform from result
        
        
        '''
        
    
    def show_waveform(self, selection=None):
        # Plot currently selected waveform
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
        rigol_control.apply_waveform(inst, 1, s, '1')

    
    def record_signals(self):
        # Record for n seconds
        inst = self.rm.open_resource(self.scope.get())
        total_time = 10
        # siglent_control.record_signals(inst, total_time)
        # V, t, ft = siglent_control.download_signals(inst, total_time)
        
        df = siglent_control.impedance(inst, total_time)
        
        # print(df)
        self.ax.clear()
        for i in df:
            # print(df[i])
            # print(df[i]['freqs'])
            self.ax.plot(df[i]['freqs'],
                         np.abs(df[i]['Z']))
            self.ax.set_xlabel('Frequency/ Hz')
            self.ax.set_xscale('log')
            self.ax.set_ylabel('Z/ $\Omega$')
            self.fig.tight_layout()
            self.canvas.draw_idle()
        
        
        


    

root = tk.Tk()
gui = MainWindow(root)
root.mainloop()