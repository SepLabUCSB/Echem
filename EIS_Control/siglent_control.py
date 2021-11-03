import pyvisa
from array import array
import numpy as np
import matplotlib.pyplot as plt
from time import sleep
import time
import pandas as pd



def get_frame_time(ftim):
    ftim = ftim.replace(' ', '')
    hr, m, s = ftim.split(':')
    t = 3600*float(hr) + 60*float(m) + float(s)
    return t   
    
    
    
def record_signals(inst, total_time):
    inst.write('TRMD STOP')
    # inst.write('MSIZ 140K')
    # inst.write('TDIV 100MS')
    start_time = time.time()
    inst.write('ARM')
    sleep(1)
    inst.write('TRMD NORMAL')
    while time.time() - start_time < total_time:
        sleep(0.001)
    inst.write('STOP')

    
 
    
def download_signals(inst, total_time, sample_time = 1):
    # Oscilloscope needs to manually be set to "history" menu
    print('Press "history" on oscilloscope')
    print('Press any key when ready')
    input('>> ')
    # print('Starting in:')
    # print('5')
    # sleep(1)
    # print('4')
    # sleep(1)
    # print('3')
    # sleep(1)
    # print('2')
    # sleep(1)
    # print('1')
    # sleep(1)
    
    
    # Determine number of frames
    inst.write('FRAM 1')
    ftim = inst.query('FTIM?')[5:-1]
    t1 = get_frame_time(ftim)
    
    inst.write('FRAM 2')
    ftim = inst.query('FTIM?')[5:-1]
    t2 = get_frame_time(ftim)
    
    frame_time = t2 - t1
    number_of_frames = int(np.floor(total_time/frame_time))
    
    
    
    
    # Store V and t as arrays, each FT to dict for each channel
    C1 = {}
    C2 = {}
    
    C1['V'] = np.array([])
    C1['t'] = np.array([])
    C1['freqs'] = {}
    C1['ft'] = {}
    
    C2['V'] = np.array([])
    C2['t'] = np.array([])
    C2['freqs'] = {}
    C2['ft'] = {}
    
    
    # Download each frame from scope, score V and t, and Fourier transform
    for frame in range(1, number_of_frames + 1):
        # Channel 1
        # Get some parameters to calculate V and t from binary
        vdiv       = float(inst.query('C1:VDIV?')[8:-2])
        voffset    = float(inst.query('C1:OFST?')[8:-2])
        sara       = float(inst.query('SARA?')[5:-5])
        
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
        C1['V'] = np.append(C1['V'], volts)
        C1['t'] = np.append(C1['t'], times)
        
        # Calculate Fourier transform. Only use first sample_time seconds
        if sample_time:
            end = np.where(times == start_time + sample_time)[0][0]
        else:
            end = None
        C1['freqs'][frame] = np.fft.rfftfreq(len(volts[:end]))*sara
        C1['ft'][frame]  = np.fft.rfft(volts[:end])
        
        
        
        # Channel 2, bad coder...
        # Get some parameters to calculate V and t from binary
        vdiv       = float(inst.query('C2:VDIV?')[8:-2])
        voffset    = float(inst.query('C2:OFST?')[8:-2])
        sara       = float(inst.query('SARA?')[5:-5])
        
        # Get starting time
        inst.write('FRAM %s'%frame)
        start_time = get_frame_time(inst.query('FTIM?')[5:-1])
        
        # Download data
        inst.write('C2:WF? DAT2')
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
        C2['V'] = np.append(C2['V'], volts)
        C2['t'] = np.append(C2['t'], times)
        
        # Calculate Fourier transform. Only use first sample_time seconds
        if sample_time:
            end = np.where(times == start_time + sample_time)[0][0]
        else:
            end = None
        C2['freqs'][frame] = np.fft.rfftfreq(len(volts[:end]))*sara
        C2['ft'][frame]  = np.fft.rfft(volts[:end])
        
        
    
        
    return C1, C2
        


def impedance(inst, total_time, sample_time=1):
    
    record_signals(inst, total_time)
    
    C1, C2 = download_signals(inst, total_time)
        
    Z = {}
    freqs = {}
    
    df = {}
    
    for i in C1['ft']:
        V = C1['ft'][i]
        I = C2['ft'][i]
        Z[i] = V/I
        freqs[i] = C1['freqs'][i]
        
        df[i] = pd.DataFrame({'Z':Z[i],
                             'freqs':freqs[i]})
        
    
    return df
    
        
        
        

    
    
    

if __name__ == '__main__':
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource('USB0::0xF4ED::0xEE3A::SDS1EDEX5R5381::INSTR')

    print(inst.query('*IDN?'))
    
    # record_signals(inst, 10)
    # C1, C2 = download_signals(inst, 10)
    #%%
    C1, C2, df = impedance(inst, 10)

