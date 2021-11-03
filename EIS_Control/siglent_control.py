import pyvisa
from array import array
import numpy as np
import matplotlib.pyplot as plt
from time import sleep
import time

rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0xF4ED::0xEE3A::SDS1EDEX5R5381::INSTR')


print(inst.query('*IDN?'))


vdiv       = float(inst.query('C1:VDIV?')[8:-2])
voffset    = float(inst.query('C1:OFST?')[8:-2])
tdiv       = float(inst.query('Time_DIV?')[5:-2])
trdl       = float(inst.query('TRDL?')[5:-2])
sara       = float(inst.query('SARA?')[5:-5])
screen_time = round(14*tdiv, 6)



def get_frame_time(ftim):
    ftim = ftim.replace(' ', '')
    hr, m, s = ftim.split(':')
    t = 3600*float(hr) + 60*float(m) + float(s)
    return t   
    
    
    
def record_signals(inst, total_time):
    inst.write('TRMD STOP')
    inst.write('MSIZ 140K')
    # inst.write('TDIV 100MS')
    sleep(0.5)
    start_time = time.time()
    inst.write('ARM')
    inst.write('TRMD NORMAL')
    while time.time() - start_time < total_time:
        sleep(0.001)
    inst.write('STOP')

    
 
    
def download_signals(inst, total_time, frame_time = None):
    # Oscilloscope needs to manually be set to "history" menu
    print('Press "history" on oscilloscope')
    print('Press any key when ready')
    input('>> ')
    
    
    # Determine number of frames
    inst.write('FRAM 1')
    ftim = inst.query('FTIM?')[5:-1]
    t1 = get_frame_time(ftim)
    
    inst.write('FRAM 2')
    ftim = inst.query('FTIM?')[5:-1]
    t2 = get_frame_time(ftim)
    
    frame_time = t2 - t1
    number_of_frames = int(np.floor(total_time/frame_time))
    
    
    # Get some parameters to calculate V and t from binary
    vdiv       = float(inst.query('C1:VDIV?')[8:-2])
    voffset    = float(inst.query('C1:OFST?')[8:-2])
    sara       = float(inst.query('SARA?')[5:-5])
    
    
    # Store V and t as arrays, each FT to dict
    V  = np.array([])
    t  = np.array([])
    ft = {}
    
    
    # Download each frame from scope, score V and t, and Fourier transform
    for frame in range(1, number_of_frames + 1):
        # Get starting time
        inst.write('FRAM %s'%frame)
        start_time = get_frame_time(inst.query('FTIM?')[5:-1])
        
        # Download data
        inst.write('C1:WF? DAT2')
        trace = inst.read_raw()
        wave = trace[22:-2]
        adc = np.array(array('b', wave))
        
        # Convert to real voltage
        volts = adc*(vdiv/25) - voffset
        
        # Create array of times
        times = np.zeros(len(volts))
        for i in range(len(volts)):
            times[i] = start_time + (1/sara)*i
        
        # Append arrays
        V = np.append(V, volts)
        t = np.append(t, times)
        
        # Calculate Fourier transform. Only use first frame_time seconds
        if frame_time:
            end = np.where(times == start_time + frame_time)[0][0]
        else:
            end = None
        freqs = np.fft.rfftfreq(len(volts[:end]))*sara
        amps  = np.fft.rfft(volts[:end])
        
        ft[frame] = np.array([freqs, amps])
        
    return V, t, ft
        



# record_signals(inst, 10)
V, t, ft = download_signals(inst, 10)
