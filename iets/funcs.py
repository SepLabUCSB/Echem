import pyvisa
import time
import numpy as np
from array import array


SR_SENS = {
    0: 100e-9,
    1: 300e-9,
    2: 1e-6,
    3: 3e-6,
    4: 10e-6,
    5: 30e-6,
    6: 100e-6,
    7: 300e-6,
    8: 1e-3,
    9: 3e-3,
    10: 10e-3,
    11: 30e-3,
    12: 100e-3,
    13: 300e-3,
    14: 1e0,
    15: 3e0,    
    }



def init_recording(inst):
        
    # Write scope params
    inst.write('TRMD AUTO')
    inst.write('TDIV 5MS')
    inst.write('MSIZ 7K')
    # inst.write('TRMD STOP')
    
    ps = {}
    # Save scope params
    ps['vdiv1']      = float(inst.query('C1:VDIV?')[8:-2])
    ps['voffset1']   = float(inst.query('C1:OFST?')[8:-2])
    ps['vdiv2' ]     = float(inst.query('C2:VDIV?')[8:-2])
    ps['voffset2']   = float(inst.query('C2:OFST?')[8:-2])
    ps['sara']       = float(inst.query('SARA?')[5:-5])
    tdiv             = float(inst.query('TDIV?')[5:-2])
    ps['frame_time'] = 14*tdiv
    
    return ps



def record_current(inst, bias, params, autolab_i_range,
                   lockin_sensitivity):
    # inst   = opened PyVisa resource
    # bias   = DC bias (str or float)
    # params = recording params
    # i_range= float
    
    ftime = params['frame_time']
    
    # Set DC offset on scope
    # inst.write(f'C1:OFST {bias}s')
    
    # Wait for 1 frame of data
    inst.write('TRMD STOP')
    inst.write('TRMD SINGLE')
    inst.query('INR?')
    inst.write('ARM')
    inst.write('FRTR')
    start_time = time.time()
    while time.time() - start_time < 10*ftime:
        inr = inst.query('INR?').strip('\n').split(' ')[1]
        if inr == '1':
            break
        time.sleep(0.001)
    inst.write('TRMD STOP')
    
    # Read data
    inst.write('C1:WF? DAT2')
    trace1   = inst.read_raw()
    wave1    = trace1[22:-2]
    adc1     = np.array(array('b', wave1))
    vdiv1    = params['vdiv1']
    voffset1 = params['voffset1']
    volts1   = adc1*(vdiv1/25) - voffset1 #scope volts
    
    # Get average
    volts1 = np.average(volts1)
    
    # Convert to lockin signal
    # Output = (signal/sensitivity - offset) x Expand x 10 V
    signal = (volts1/10)*lockin_sensitivity
    
    # Convert to autolab current
    i = signal*autolab_i_range
    
    return i



def record_current_srs(srs, autolab_i_range, recording_time):
    
    st    = time.time()
    vs    = []    

    # srs.write('AGAN')    
        
    while (time.time() - st) < recording_time:
        
        v = float(srs.query('OUTP? 3').strip('\n'))
        vs.append(v)
        
    v = np.average(vs)
    i = v*autolab_i_range
    return i


def record_current_sr(sr, sr_sens, autolab_i_range, recording_time):
    
    st    = time.time()
    vs    = []    

    # srs.write('AGAN')    
        
    while (time.time() - st) < recording_time:
        
        v = float(sr.query('MAG').strip('\n').strip('\r'))
        vs.append(v)
    
        
    v = np.average(vs)
    sens = SR_SENS[sr_sens]
    v = sens*v/10000
    i = v*autolab_i_range
    
    return i


if __name__ == '__main__':
    
    # scope   = 'USB0::0xF4ED::0xEE3A::SDS1EDED5R0471::INSTR'
    scope = 'GPIB0::16::INSTR'
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(scope)
    
    
    # initialize scope params
    # params = init_recording(inst)
    
    ts = []
    Is = []
    for _ in range(100):
        st = time.time()
        i = record_current_sr(inst, 9, 1, 0.5)
        t = time.time()-st
        
        print(t, i)
        
        ts.append(t)
        Is.append(i)
        
    print('')
    print('5 ms tdiv:')
    print(f'Average time: {np.average(ts):.05f} +- {np.std(ts):.05f}')
    print(f'Average I: {np.average(Is):.05f} +- {np.std(Is):.05f}')
    

