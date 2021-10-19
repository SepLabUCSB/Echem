import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

sample_freq = 100000
total_time = 10

def create_waveform(freqs, phases, sample_freq,
                    total_time, amax=0.01, amps = None,
                    save_path = None, plot = True):
    
    time = np.arange(0, total_time, 1/sample_freq)
    
    S = 0   # Summed multisin waveform    
    i = 0   
    
    
    amax = 2*amax #correct HEKA 0.5*voltage factor
    

    if amps is None:
        amps = np.ones(len(freqs))
    

    while i < len(freqs): #add up sines with freq, phase, amplitude from csv file
        sin = amps[i] * np.sin(2 * np.pi * freqs[i] * time + phases[i])
        S = S + sin
        i += 1
            
            
    S = S * (amax/S.max()) #normalize maximum potential to amax
    S = np.array(S)
    
    if save_path:
        file = os.path.join(save_path, 'waveform.asc')
        pd.Series(S).to_csv(file, header=None, index=False,
                            encoding='ascii')
    
    else:
        file = os.path.join(r'C:\Users\BRoehrich\Desktop', 'waveform.asc')
        pd.Series(S).to_csv(file, header=None, index=False,
                            encoding='ascii')
    
    
    if plot:
        fig, ax = plt.subplots()
        ax.plot(time, S)
        ax.set_xlabel('Time/ s')
        ax.set_ylabel('Amplitude/ V')
        
        fig, ax = plt.subplots()
        ax.plot(sample_freq*np.fft.rfftfreq(len(S))[0:],
                abs(np.fft.rfft(S)[0:]))
        ax.set_xlabel('Frequency/ Hz')
        ax.set_ylabel('Amplitude')
        ax.set_xscale('log')
    
    print('')
    print('Created new waveform')
    print('')
    print('Max frequency: ', max(freqs))
    print('Minimum frequency: ', min(freqs))
    print('Amplitude: ', amax/2)
    print('Sampling frequency: ', sample_freq)
    print('Total time: ', total_time)
    print('Saved to ', os.path.join(r'C:\Users\BRoehrich\Desktop', 'waveform.asc'))
    
    
    return S


def waveform_from_result(FT_EIS, sample_freq, total_time, 
                         amax=0.01, save_path = None, plot = True):
    '''
    Generates new waveform with amplitudes determined from
    
    given EIS spectrum. A_j proportional to |Z_j|
    
    '''    
    df = FT_EIS.ft[1]
    
    freqs = df['f'].to_numpy()
    amps = np.sqrt(np.absolute(df['Z'].to_numpy()))
    
    # Find file with frequency array
    file = FT_EIS.csv_dir + r'\f_%s_%s.csv' %(
        FT_EIS.highest_freq, FT_EIS.lowest_freq)
    
    d = pd.read_csv(file, skiprows=1, names=('index', 'freqs', 'phases'),
                    dtype=float)
    file_freqs = d['freqs'].to_numpy()
    phases = d['phases'].to_numpy()
    
    
    
    if file_freqs.all() != freqs.all():
        print("ERROR: Frequencies don't match!")
        print('Frequencies from Fourier transform:')
        print(d['freqs'].to_numpy())
        print('')
        print('Frequencies from csv file:')
        print(freqs)
        
        return None
    
    
    else:
        S = create_waveform(freqs, phases, sample_freq, total_time, amax,
                        amps=amps, save_path=save_path, plot=plot)
        
        
        return S
        
    
    
    
    
    
    
    
