import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import numpy as np
from scipy import signal
import pyvisa
import time
from time import sleep
from array import array
from scipy.signal import butter, lfilter, freqz
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


class MainWindow:
    
    global this_dir, rigol_waves
    
    def __init__(self, root):
        self.root = root
        root.title("Fast SWV")
        root.attributes('-topmost', 1)
        root.attributes('-topmost', 0)    

        self.rm = pyvisa.ResourceManager()
            
        
        # Initialize frames and canvas
        
        # frame:  left
        self.frame = tk.Frame(self.root)
        self.frame.grid(row=0, column = 0)
        
        # fig: right
        self.fig = plt.Figure(figsize=(5,4), dpi=100)
        self.ax  = self.fig.add_subplot(111)
                
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().grid(row=0, column=1)
        
        
        
        ########################################
        ###         SWV Parameters           ###
        ########################################
        
        # E initial
        text = tk.Label(self.frame, text='E_initial (V):')
        text.grid(row=0, column = 1, sticky = 'E')
        self.E_initial = tk.Text(self.frame, height=1, width=7)
        self.E_initial.insert('1.0', '0.0')
        self.E_initial.grid(row=0, column=2)
        
        # E final
        text = tk.Label(self.frame, text='E_final (V):')
        text.grid(row=1, column = 1, sticky = 'E')
        self.E_final = tk.Text(self.frame, height=1, width=7)
        self.E_final.insert('1.0', '-0.2')
        self.E_final.grid(row=1, column=2)
        
        # Pulse height
        text = tk.Label(self.frame, text='Pulse height (V):')
        text.grid(row=2, column = 1, sticky = 'E')
        self.P_height = tk.Text(self.frame, height=1, width=7)
        self.P_height.insert('1.0', '0.025')
        self.P_height.grid(row=2, column=2)
        
        # Step height
        text = tk.Label(self.frame, text='Step height (V):')
        text.grid(row=3, column = 1, sticky = 'E')
        self.S_height = tk.Text(self.frame, height=1, width=7)
        self.S_height.insert('1.0', '0.005')
        self.S_height.grid(row=3, column=2)
        
        # Pulse width
        text = tk.Label(self.frame, text='Pulse width (s):')
        text.grid(row=4, column = 1, sticky = 'E')
        self.P_width = tk.Text(self.frame, height=1, width=7)
        self.P_width.insert('1.0', '0.1')
        self.P_width.grid(row=4, column=2)
        
        # Current range
        text = tk.Label(self.frame, text='Current range:')
        text.grid(row=5, column = 1, sticky = 'E')
        self.current_range = tk.Text(self.frame, height=1, width=7)
        self.current_range.insert('1.0', '1e-6')
        self.current_range.grid(row=5, column=2)

        

        # AWG Control
        text = tk.Label(self.frame, text='Arb:')
        text.grid(row=6, column=0)
        self.arb = tk.StringVar(self.frame)
        
        try:
            # Look for Rigol arb, i.e.
            # 'USB0::0x1AB1::0x0643::DG8A232302748::INSTR'
            default_arb = [inst for inst in self.rm.list_resources() 
                            if len(inst.split('::')) > 3 
                            and inst.split('::')[3].startswith('DG')][0]
        
        except:
            default_arb = ''
        
        self.arb.set(default_arb)
        self.arb_selector = tk.OptionMenu(self.frame, self.arb, 
                                               *self.rm.list_resources())
        self.arb_selector.grid(row=6, column=1)
        self.apply_waveform_button = tk.Button(self.frame, text='Apply Wave', 
                                               command=self.apply_waveform)
        self.apply_waveform_button.grid(row=6, column=2)
        
        
        
        # Scope control
        
        text = tk.Label(self.frame, text='Scope:')
        text.grid(row=7, column=0)
        self.scope = tk.StringVar(self.frame)
        
        try:
            default_scope = [inst for inst in self.rm.list_resources() 
                            if len(inst.split('::')) > 3 
                            and inst.split('::')[3].startswith('SDS')][0]
        
        except:
            default_scope = ''
        
        self.scope.set(default_scope)
        self.scope_selector = tk.OptionMenu(self.frame, self.scope, 
                                               *self.rm.list_resources())
        self.scope_selector.grid(row=7, column=1)
        
        self.record_signals_button = tk.Button(self.frame, text='Record', 
                                               command=self.record_frame)
        self.record_signals_button.grid(row=7, column=2)
        
    
    
    
    def apply_waveform(self):
        '''
        Apply arbitrary waveform s to Rigol DG812 AWG
        
        https://rigol.force.com/support/s/article/methods-for-programmatically-creating-arbitrary-waves1
    
        Parameters
        ----------
        inst : pyvisa Instrument
            Handle for DG812.
            USB0::0x1AB1::0x0643::DG8A232302748::INSTR
        channel: 1 or 2
            Channel to apply waveform to
        s : array of int16
            Arbitrary waveform to apply.
        '''
        
        # Create staircase waveform from user parameters
        
        E_1 = float(self.E_initial.get('1.0', 'end'))
        E_2 = float(self.E_final.get('1.0', 'end'))
        p_h = float(self.P_height.get('1.0', 'end'))
        s_h = float(self.S_height.get('1.0', 'end'))
        p_w = float(self.P_width.get('1.0', 'end'))
        freq = 5000
                
        if E_1 > E_2:
            p_h = -p_h
            s_h = -s_h
    
        # Every 2*p_w seconds, gain s_h volts
        # Total time is 2*p_w*|E_2 - E_1|/s_h
        
        # print(p_w, E_1, E_2, s_h)
        
        total_time = 2*p_w* np.abs(E_2-E_1) / np.abs(s_h)
        t = np.arange(0, total_time, 1/freq)
        # print(t)
        
        # Initialize V
        V = np.zeros(len(t))
        
        # Determine staircase voltammetry values
        steps = np.arange(E_1, E_2, s_h)
        self.number_of_steps = len(steps)
        
        # Find # of points per staircase step
        pts_per_step = int(len(V)/len(steps))
        
        
        # Staircase waveform
        for i in range(len(steps)):
            V[pts_per_step*i : pts_per_step*(i+1)] = steps[i]
        
        sw = p_h*signal.square(2*np.pi*(1/(2*p_w))*t)/2
                
        s = sw + V
        
        # Starting and ending hold potentials
        start = E_1*np.ones(3*freq)
        end = E_2*np.ones(3*freq)
        
        s = np.concatenate((start, s, end))
        
        
        self.frame_time = len(s)/freq
    
        
        # Send waveform to arb
        channel = 1
        Vpp = 4*(max(abs(s)) - E_1)
        s = to_int16(s)
        
        
        inst = self.rm.open_resource(self.arb.get())
        
        try:
            inst.write('*RST')
            inst.timeout = 1000
            inst.write(':SOURCE%s:APPL:SEQ'%(channel))
            wait(inst)
            inst.write(':SOURCE%s:FUNC:SEQ:FILT INSERT'%(channel))
            wait(inst)
            blocks = send_bytes(inst, s, channel)
            wait(inst)
            inst.write(':SOURCE%s:VOLTAGE %sVPP'%(channel, Vpp))
            inst.write(':SOURCE%s:VOLTAGE:OFFSET 0'%(channel))
            inst.write(':SOURCE%s:FUNC:SEQ:SRAT %s'%(channel, freq))
            inst.write(':SOURCE%s:FUNC:SEQ:EDGETime 0.000005'%(channel))
            inst.write(':SOURCE1:FUNC:SEQ:PER 1 0')
            inst.write(':SOURCE%s:FUNC:'%(channel))
            wait(inst)
            # inst.write(':OUTPUT%s ON;'%(channel))
            # wait(inst)
            print('Waveform loaded!\n')
            inst.clear()
        
        except:
            print('Could not connect')
            pass
        
        return blocks
    
    
    def record_frame(self):
        inst = self.rm.open_resource(self.scope.get())
        arb  = self.rm.open_resource(self.arb.get())
        
        inst.write('TRMD AUTO') 
        
        frame_time = self.frame_time
        
        if frame_time < 7:
            inst.write('TDIV 1s')
        elif frame_time >= 7 and frame_time < 14:
            inst.write('TDIV 1s')
        elif frame_time >= 14 and frame_time < 28:
            inst.write('TDIV 2s')
        elif frame_time >= 28 and frame_time < 70:
            inst.write('TDIV 5s')
        elif frame_time >= 70 and frame_time < 140:
            inst.write('TDIV 10s')
        elif frame_time >= 140 and frame_time < 280:
            inst.write('TDIV 20s')
        elif frame_time >= 280 and frame_time < 700:
            inst.write('TDIV 50s')
        elif frame_time >= 700:
            inst.write('TDIV 100s')
        
        
        
        vdiv1       = float(inst.query('C1:VDIV?')[8:-2])
        voffset1    = float(inst.query('C1:OFST?')[8:-2])
        vdiv2       = float(inst.query('C2:VDIV?')[8:-2])
        voffset2    = float(inst.query('C2:OFST?')[8:-2])
        tdiv        = float(inst.query('TDIV?')[5:-2])
        self.sara   = float(inst.query('SARA?')[5:-5])
        
        if self.sara < 5000:
            print('Increase memory depth! Sampling rate too low.')
            return 0,0
        elif self.sara > 50000:
            print('Decrease memory depth! Sampling rate too high.')
            return 0,0

        
        # Turn on AWG and start recording
        arb.write(':OUTPUT1 OFF')
        # arb.write(':OUTPUT2 OFF')
        wait(arb)
        
        arb.write(':OUTPUT1 ON')   
        # arb.write(':OUTPUT2 ON')
        
        frame_start_time = time.time()
        
        inst.write('TRMD SINGLE;ARM;FRTR ')
        print('Starting')
        
        while time.time() - frame_start_time < frame_time-1:
            time.sleep(0.01) 
            
        inst.write('STOP')
        arb.write(':OUTPUT1 OFF')
        print(f'Aquisition complete, {time.time() - frame_start_time} s')
            
        inst.timeout = 10000
        # Get CH 1 data
        inst.write('C1:WF? DAT2')
        # wait(inst)
        trace1 = inst.read_raw()
        wave1 = trace1[22:-2]
        adc1 = np.array(array('b', wave1))
        
        # Get CH 2 data
        inst.write('C2:WF? DAT2')
        # wait(inst)
        trace2 = inst.read_raw()
        wave2 = trace2[22:-2]
        adc2 = np.array(array('b', wave2))
        
        # Convert to voltages
        volts1 = adc1*(vdiv1/25) - voffset1 
        volts2 = adc2*(vdiv2/25) - voffset2 
        
        self.volts1 = volts1
        self.volts2 = volts2
        
        self.process_data(volts1, volts2)
        
        
        
        
    
    def process_data(self, volts1, volts2):
        # volts1: CH1 data (Voltage)
        # volts2: CH2 data (Current)
        
        E = -volts1
        I = volts2/float(self.current_range.get('1.0'))
        
        E = butter_lowpass_filter(E, 2000, self.sara, order=5)
        I = butter_lowpass_filter(I, 2000, self.sara, order=5)
        
        E = E[100:]
        I = I[100:]
               
        # fig, ax = plt.subplots(figsize=(3,3), dpi=150)
        # ax.plot(I)
        # ax.plot(np.abs(np.gradient(I))) 
        
        p_w = float(self.P_width.get('1.0', 'end'))
        
        # Figure out where the staircase waveform starts in the sketchiest way possible
        avg_grad_I = max(abs((np.gradient(I[:int(self.sara)]))))
        print(avg_grad_I)
        start_index = np.where(np.abs(np.gradient(I)) > 
                               4*np.abs(avg_grad_I))[0][0]
        start_index = int(start_index - p_w*self.sara)
        
        
        
        E = E[start_index:]
        I = I[start_index:]
             
        
        
        
                
        
        steps_v = []
        steps_i = []
        
        
        for j in range(2*self.number_of_steps+2):
            # Break up into each square wave step
            a = int(self.sara*p_w*j)
            b = int(self.sara*p_w*(j+1))
            
            
            steps_v.append(E[a:b])
            steps_i.append(I[a:b])
        
        fig, ax = plt.subplots(figsize=(3,3), dpi=150)
        plt.plot(np.concatenate([arr for arr in steps_v]))
        self.ax.clear()
                
        for sw_freq in [750,500,100,10]:
            I_delta, V_delta, V_step = square_wave(steps_v, steps_i, sw_freq,
                                                   self.sara, p_w)
            
            x = V_step[1::2]
            y = I_delta[1::2]
            
            # x = V_step
            # y = I_delta/sw_freq
            
            if sw_freq == 10:
                self.V_step = V_step
                self.I_delta = I_delta
                self.V_delta = V_delta
            
            self.ax.plot(x, y, label=f'{sw_freq}')
        
        
        self.ax.set_xlabel('$V_{step}$ / V vs SCE')
        self.ax.set_ylabel('$\delta$I / $\mu$A')
        self.ax.legend()
        self.fig.tight_layout()
        self.canvas.draw_idle()
        
        freqs = self.sara*np.fft.rfftfreq(len(steps_v[0]))[1:]
        freqs = freqs[freqs <= 1000]
        sw_i = np.array([])
        
        for sw_freq in freqs:
            I_delta, V_delta, self.V_step = square_wave(steps_v, steps_i, sw_freq,
                                                   self.sara, p_w)
            if len(sw_i) == 0:
                sw_i = I_delta[::2]
            else:
                sw_i = np.vstack((sw_i, I_delta[::2]))
        
        self.freqs = freqs
        self.sw_i = sw_i.T
        self.V_step = V_step[::2] # Select forward steps
        
        fig, ax = plt.subplots(figsize=(3,3), dpi=150)
        plot = plt.pcolormesh(freqs, self.V_step, -self.sw_i/freqs, cmap='RdBu_r', shading='gouraud')
        plt.xscale('log')
        cbar = plt.colorbar(plot)
        cbar.set_label('$\delta$I/f /mC')
        plt.xlabel('Frequency/ Hz')
        plt.ylabel('Potential/ V vs SCE')
        
    
        
        


# Data transfer functions
def wait(inst):
    r = False
    
    while r is False:
        try:
            r = inst.query('*OPC?')
        except:
            sleep(0.1)

            
def to_int16(signal):
    signal = np.array(signal)
    
    if signal.max() > abs(signal.min()):
    
        # signal = signal - signal.min()
        signal = np.int16(((signal/signal.max()) * int("3fff", 16)))
    
    elif abs(signal.min()) > signal.max():
        
        signal = np.int16(((signal/abs(signal.min())) * int("3fff", 16)))
        
    elif abs(signal.min()) == signal.max():
        signal = np.int16(((signal/signal.max()) * int("3fff", 16)))
    
    return signal


def send_bytes(inst, signal, channel):
    # Break signal up into 16kpts chunks (32 kB)
    number_of_blocks = int(np.floor(len(signal)/16000))
    blocks = {}
    
    string = ':SOURCE%s:TRACe:DATA:DAC16 VOLATILE,CON,'%channel
    end_string = ':SOURCE%s:TRACe:DATA:DAC16 VOLATILE,END,'%channel
    
    
    if len(signal)%16000 == 0:
        number_of_blocks -= 1
    
    for i in range(number_of_blocks+1):
        blocks[i] = signal[16000*i:16000*(i+1)]
    
    for i in range(number_of_blocks):
        # print('Sending points %s:%s'%(blocks[i][0], blocks[i][-1]))
        inst.write_binary_values(string, blocks[i], datatype='h')
        wait(inst)
    
    inst.write_binary_values(end_string, blocks[number_of_blocks], datatype='h')
    # print('Sending points %s:%s'%(blocks[number_of_blocks][0], blocks[number_of_blocks][-1]))
    
    return blocks


def square_wave(V, I, sw_freq, sample_freq, step_duration):
    
    # Get index to sample current from
    # Assumes each step has a total duration of 1 s
    ind1 = int(0.9*sample_freq / sw_freq)
    ind2 = int(sample_freq / sw_freq)

    # print(ind1, ind2) 
    
    I_delta = []
    V_delta = []
    V_step  = []
    
    for j in range(1, len(I)):
        try:
            delta_I = np.average(I[j][ind1:ind2] - I[j-1][ind1:ind2])
        except:
            print(j, I[j], I[j-1])
            print(len(I), sw_freq)
        delta_V = np.average(V[j][ind1:ind2] - V[j-1][ind1:ind2])

        I_delta.append(delta_I)
        V_delta.append(delta_V)
        
        V_step.append(np.average(V[j][ind1:ind2]))
    
    return np.array(I_delta), np.array(V_delta), np.array(V_step)






root = tk.Tk()
gui = MainWindow(root)
root.mainloop()



