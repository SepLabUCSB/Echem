import pyvisa
from array import array
import numpy as np
import matplotlib.pyplot as plt
from time import sleep
import time
import pandas as pd


class FourierTransformData:
    
    def __init__(self, time, freqs, CH1data, CH2data, Z = None,
                 phase = None, waveform = None):
        self.time = time
        self.freqs = freqs
        self.CH1data = CH1data
        self.CH2data = CH2data
        
        self.Z = Z
        self.phase = phase
        self.waveform = waveform




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
    
        

def record_continuous(inst, total_time, sample_time=1):
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
    
    
    V1      = np.array([])
    V2      = np.array([])
    t       = np.array([])
    ft      = {}
    
    frame_time = 14*tdiv
    start_time = time.time()
       
    frame = 0
    
    while time.time() - start_time < total_time:
        
        # Determine t=0 for frame
        frame_start_time = time.time()
       
        # Record frame
        inst.write('TRMD AUTO')
        sleep(1.2*frame_time)
        
        # Get CH 1 data
        inst.write('C1:WF? DAT2')
        trace1 = inst.read_raw()
        wave1 = trace1[22:-2]
        adc1 = np.array(array('b', wave1))
        
        # Get CH 2 data
        inst.write('C2:WF? DAT2')
        trace2 = inst.read_raw()
        wave2 = trace2[22:-2]
        adc2 = np.array(array('b', wave2))
        
        # Convert to voltages
        volts1 = adc1*(vdiv1/25) - voffset1 
        volts2 = adc2*(vdiv2/25) - voffset2  
        
        # Get time array
        times = np.zeros(len(volts1))
        for i in range(len(volts1)):
            times[i] = frame_start_time + (1/sara)*i - start_time
           
        # Save results
        V1 = np.append(V1, volts1)
        V2 = np.append(V2, volts2)
        t  = np.append(t, times)
        
        # Only Fourier transform first sample_time s
        if sample_time:
            end = np.where(times == times[0] + sample_time)[0][0]
        else:
            end = None
        
        freqs = sara*np.fft.rfftfreq(len(volts1[:end]))[1:]
        ft1   =      np.fft.rfft(volts1[:end])[1:]
        ft2   =      np.fft.rfft(volts2[:end])[1:]
        
        ft[frame] = FourierTransformData(time    = times[0],
                                         freqs   = freqs,
                                         CH1data = ft1,
                                         CH2data = ft2,)
        
        
        frame += 1
        
    
    return ft
        


def record_single(inst, start_time, frame_time,
                  vdiv1, voffset1, vdiv2, voffset2,
                  sara, sample_time=1):
    # Determine t=0 for frame
    frame_start_time = time.time()
   
    # Record frame
    inst.write('TRMD AUTO')
    sleep(1.2*frame_time)
    
    # Get CH 1 data
    inst.write('C1:WF? DAT2')
    trace1 = inst.read_raw()
    wave1 = trace1[22:-2]
    adc1 = np.array(array('b', wave1))
    
    # Get CH 2 data
    inst.write('C2:WF? DAT2')
    trace2 = inst.read_raw()
    wave2 = trace2[22:-2]
    adc2 = np.array(array('b', wave2))
    
    # Convert to voltages
    volts1 = adc1*(vdiv1/25) - voffset1 
    volts2 = adc2*(vdiv2/25) - voffset2  
    
    # Get time array
    times = np.zeros(len(volts1))
    for i in range(len(volts1)):
        times[i] = frame_start_time + (1/sara)*i - start_time
           
    # Only Fourier transform first sample_time s
    if sample_time:
        end = np.where(times == times[0] + sample_time)[0][0]
    else:
        end = None
    
    freqs = sara*np.fft.rfftfreq(len(volts1[:end]))[1:]
    ft1   =      np.fft.rfft(volts1[:end])[1:]
    ft2   =      np.fft.rfft(volts2[:end])[1:]
    
    ft = FourierTransformData(time    = times[0],
                              freqs   = freqs,
                              CH1data = ft1,
                              CH2data = ft2,)
    
    return ft

    
    

if __name__ == '__main__':
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource('USB0::0xF4ED::0xEE3A::SDS1EDEX5R5381::INSTR')

    print(inst.query('*IDN?'))
    
    ft = record_continuous(inst, 10)
    

