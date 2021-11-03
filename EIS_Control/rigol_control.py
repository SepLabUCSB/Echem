import matplotlib.pyplot as plt
import numpy as np
import pyvisa
from time import sleep
import pandas as pd
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')



def to_int16(signal):
    
    signal = np.array(signal)
    
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
    


def send_bytes(inst, signal, channel):
    # Break signal up into 16kpts chunks (32 kB)
    number_of_blocks = int(np.floor(len(signal)/16000))
    blocks = {}
    
    string = ':SOURCE%s:TRACe:DATA:DAC16 VOLATILE,CON,'%channel
    end_string = ':SOURCE%s:TRACe:DATA:DAC16 VOLATILE,END,'%channel
    
    for i in range(number_of_blocks+1):
        blocks[i] = signal[16000*i:16000*(i+1)]

    
    for i in range(number_of_blocks):
        # print('Sending points %s:%s'%(blocks[i][0], blocks[i][-1]))
        inst.write_binary_values(string, blocks[i], datatype='h')
        wait(inst)
    
    inst.write_binary_values(end_string, blocks[number_of_blocks], datatype='h')
    # print('Sending points %s:%s'%(blocks[number_of_blocks][0], blocks[number_of_blocks][-1]))
    
    return blocks



def wait(inst):
    r = False
    
    while r is False:
        try:
            r = inst.query('*OPC?')
        except:
            sleep(0.2)
    


def apply_waveform(inst, channel, s, Vpp):
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
    
    if not s.dtype == 'int16':
        s = to_int16(s)
    
    
    if channel == 1 or channel == 2:
        channel = int(channel)
    
    else:
        print('Invalid channel. Must be 1 or 2.')
        import sys
        sys.exit()
    
    
    try:
        print(inst.query('*IDN?'))
    except:
        print('Could not connect')
        import sys
        sys.exit()
    
    # inst.write('*RST')
    inst.write(':SOURCE%s:APPL:SEQ'%(channel))
    wait(inst)
    inst.write(':SOURCE%s:FUNC:SEQ:FILT INSERT'%(channel))
    wait(inst)
    send_bytes(inst, s, channel)
    wait(inst)
    inst.write(':SOURCE%s:VOLTAGE %sVPP'%(channel, Vpp))
    inst.write(':SOURCE%s:FUNC:SEQ:SRAT 100000'%(channel))
    inst.write(':SOURCE%s:FUNC:SEQ:EDGETime 0.000005'%(channel))
    # inst.write(':SOURCE1:FUNC:SEQ:PER 1 1')
    inst.write(':SOURCE%s:FUNC:'%(channel))
    wait(inst)
    inst.write(':OUTPUT%s ON;'%(channel))
    wait(inst)
    print('Waveform loaded!')
    inst.clear()


#%%

if __name__ == '__main__':
    Npts = 100000

    # t = np.arange(0,Npts)
    
    # w1 = 0.5/100000
    # w2 = 1/100000
    # s = np.sin(2*np.pi*t*w1) + np.sin(2*np.pi*t*w2 + np.pi/4)
    
    df = pd.read_csv(r'C:/Users/BRoehrich/Desktop/Waveforms/Rigol_100k_1k_1_16freqs.csv', 
                      skiprows=9, names=('x', 'V'))
    
    s = df['V'].to_numpy()
    
    
    rm = pyvisa.ResourceManager()
    rigol = rm.open_resource('USB0::0x1AB1::0x0643::DG8A232302748::INSTR')
    
    rigol.write('*RST')
    apply_waveform(rigol, 1, to_int16(s), '1')
    apply_waveform(rigol, 2, to_int16(s), '1')

    
        




