import pyvisa
import time
import numpy as np
from array import array



def init_recording(scope):
    
    # Connect to scope
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(scope)
    
    # Write scope params
    inst.write('TRMD AUTO')
    inst.write('TDIV 50MS')
    inst.write('MSIZ 7K')
    inst.write('TRMD STOP')
    
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



def record_current(inst, bias, params):
    # inst   = opened PyVisa resource
    # bias   = DC bias (str or float)
    # params = recording params
    
    ftime = params['frame_time']
    
    # Set DC offset on scope
    inst.write(f'C1:OFST {bias}s')
    
    # Wait for 1 frame of data
    inst.write('TRMD STOP')
    inst.write('ARM')
    inst.write('FRTR')
    start_time = time.time()
    while time.time() - start_time < 1.2*ftime:
        time.sleep(0.001)
    inst.write('TRMD STOP')
    
    # Read data
    inst.write('C1:WF? DAT2')
    trace1   = inst.read_raw()
    wave1    = trace1[22:-2]
    adc1     = np.array(array('b', wave1))
    vdiv1    = params['vdiv1']
    voffset1 = params['voffset1']
    volts1   = adc1*(vdiv1/25) - voffset1
    
    # Divide by current range
    
    # i = avg()
    
    return #i



