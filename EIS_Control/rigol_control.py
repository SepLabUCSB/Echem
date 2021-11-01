import matplotlib.pyplot as plt
import numpy as np
import pyvisa
from time import sleep
import pandas as pd
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')
fileout = "gen.RAF"



# Npts = 100000

# t = np.arange(0,Npts)

# w1 = 10/Npts
# w2 = 6/Npts
# s = np.sin(2*np.pi*t*w1) + np.sin(2*np.pi*t*w2 + np.pi/4)

df = pd.read_csv(r'C:/Users/BRoehrich/Desktop/Waveforms/Rigol_100k_10k_1_44freqs.csv', 
                 skiprows=9, names=('x', 'V'))

s = df['V'].to_numpy()


def to_int16(signal):
    
    if signal.max() > abs(signal.min()):
    
        # signal = signal - signal.min()
        signal = np.int16(((signal/signal.max()) * int("3fff", 16)))
    
    elif abs(signal.min()) > signal.max():
        
        signal = np.int16(((signal/abs(signal.min())) * int("3fff", 16)))
    
    
    return signal



def get_bytes(val):
    # Unused but saving just in case...
    bit = int(8191.5*val + 8191.5)
    
    bigbyte = int(bit/2**8)
    littlebyte = bit-(bigbyte*2**8)
    return '%c%c'%(littlebyte, bigbyte)
    


def send_bytes(inst, signal):
    # Break signal up into 16kpts chunks (32 kB)
    number_of_blocks = int(np.floor(len(signal)/16000))
    blocks = {}
    
    string = ':SOURCE1:TRACe:DATA:DAC16 VOLATILE,CON,'
    end_string = ':SOURCE1:TRACe:DATA:DAC16 VOLATILE,END,'
    
    for i in range(number_of_blocks+1):
        blocks[i] = signal[16000*i:16000*(i+1)]

    
    for i in range(number_of_blocks):
        # print('Sending points %s:%s'%(blocks[i][0], blocks[i][-1]))
        inst.write_binary_values(string, blocks[i], datatype='h')
        wait()
    
    inst.write_binary_values(end_string, blocks[number_of_blocks], datatype='h')
    # print('Sending points %s:%s'%(blocks[number_of_blocks][0], blocks[number_of_blocks][-1]))
    
    return blocks



def wait():
    r = False
    
    while r is False:
        try:
            r = inst.query('*OPC?')
        except:
            sleep(0.2)
    

rm = pyvisa.ResourceManager()
inst = rm.open_resource('USB0::0x1AB1::0x0643::DG8A232302748::INSTR')

def apply_waveform(inst, s):
    '''
    Apply arbitrary waveform s to Rigol DG812 AWG
    
    https://rigol.force.com/support/s/article/methods-for-programmatically-creating-arbitrary-waves1

    Parameters
    ----------
    inst : pyvisa Instrument
        Handle for DG812.
        USB0::0x1AB1::0x0643::DG8A232302748::INSTR
    s : array of int16
        Arbitrary waveform to apply.
    '''
    
    
    try:
        print(inst.query('*IDN?'))
    except:
        print('Could not connect')
        import sys
        sys.exit()
    
    inst.write('*RST')
    wait()
    inst.write(':SOURCE1:APPL:SEQ')
    wait()
    inst.write(':SOURCE1:FUNC:SEQ:FILT INSERT')
    wait()
    send_bytes(inst, s)
    wait()
    inst.write(':SOURCE1:VOLTAGE 0.02VPP')
    inst.write(':SOURCE1:FUNC:SEQ:SRAT 100000')
    inst.write(':SOURCE1:FUNC:SEQ:EDGETime 0.000005')
    # inst.write(':SOURCE1:FUNC:SEQ:PER 1 1')
    inst.write(':SOURCE1:FUNC:')
    wait()
    inst.write(':OUTPUT2 OFF')
    inst.write(':OUTPUT1 ON;')
    wait()
    inst.clear()


#%%

apply_waveform(inst, to_int16(s))
    
        




