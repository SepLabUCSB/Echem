import pyvisa
import time
import numpy as np
from array import array



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


    



if __name__ == '__main__':
    
    scope   = 'USB0::0xF4ED::0xEE3A::SDS1EDED5R0471::INSTR'
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(scope)
    
    
    # initialize scope params
    params = init_recording(inst)
    
    ts = []
    Is = []
    for _ in range(1000):
        st = time.time()
        i = record_current(inst, 0, params, 1, 10e-3)
        t = time.time()-st
        
        print(t, i)
        
        ts.append(t)
        Is.append(i)
        
    print('')
    print('5 ms tdiv:')
    print(f'Average time: {np.average(ts):.05f} +- {np.std(ts):.05f}')
    print(f'Average I: {np.average(Is):.05f} +- {np.std(Is):.05f}')
    

