import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

default_save_path = r'C:\Users\BRoehrich\Desktop'


def get_units(num):        
        if num < 1000:
            return ('', 1)
        
        if num >= 1e3:
            return('k', 1e3)



def create_waveform(freqs, phases, sample_freq,
                    total_time, amax, 
                    amps = None, csv = False, Rigol = False,
                    save_path = None, save = True, plot = True):
    
    time = np.arange(0, total_time, 1/sample_freq)
    
    freqs = np.asarray(freqs)
    phases = np.asarray(phases)
    
    S = 0   # Summed multisin waveform    
    i = 0   
    
    if not Rigol:
        amax = 2*amax #correct HEKA 0.5*voltage factor
    

    if amps is None:
        amps = np.ones(len(freqs))
    

    while i < len(freqs): #add up sines with freq, phase, amplitude from csv file
        sin = amps[i] * np.sin(2 * np.pi * freqs[i] * time + phases[i])
        S = S + sin
        i += 1
            
    
    if abs(S.max()) > abs(S.min()):
        S = S * (amax/S.max()) #normalize maximum potential to amax
    elif abs(S.min()) > abs(S.max()):
        S = S* (amax/abs(S.min()))
    S = np.array(S)
    
    
    
    
    # Save file name convention
    unit1, factor1 = get_units(sample_freq)
    unit2, factor2 = get_units(max(freqs))
    unit3, factor3 = get_units(min(freqs))
    save_name = 'EIS_%s%s_%s%s_%s%s_%sfreqs'%(round(sample_freq/factor1), unit1,
                                      round( max(freqs)/factor2), unit2, 
                                      round( min(freqs)/factor3), unit3,
                                      len(freqs))
    
    if not save_path:
        save_path = default_save_path
    
    
    if csv:
        if save:
            file = os.path.join(save_path, save_name + '.csv')
            pd.Series(S).to_csv(file, header=None, index=False)
        
            
    else:
        if save:
            file = os.path.join(save_path, save_name + '.asc')
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
        ax.set_title(save_name)
    
    
    print('')
    print('')
    print('Created new waveform')
    print('')
    print('Max frequency: %s Hz' %max(freqs))
    print('Min frequency: %s Hz' %min(freqs))
    print('Amplitude: %s mV' %(1000*S.max()))
    print('Sampling frequency: %s Hz' %sample_freq)
    print('Total time: %s s' %total_time)
    if save:
        print('Saved to ', os.path.join(save_path, save_name))
    
    
    return S, save_name




def waveform_from_result(FT_EIS, sample_freq, total_time, 
                         amax=0.01, Rigol = False, save = True,
                         save_path = None, plot = True):
    '''
    Generates new waveform with amplitudes determined from
    
    given EIS spectrum. A_j proportional to |Z_j|
    
    '''    
    df = FT_EIS.ft[1]
    
    freqs = df['f'].to_numpy()
    amps = np.sqrt(np.absolute(df['Z'].to_numpy()))
    
    # Find file with frequency array
    try:
        file = FT_EIS.csv_dir + r'\f_%s_%s_%sfreqs.csv' %(
            FT_EIS.highest_freq, FT_EIS.lowest_freq, FT_EIS.number_of_freqs)
    except:
        print('File not recognized:')
        print(os.path.join(FT_EIS.csv_dir, r'\f_%s_%s_%sfreqs.csv' %(
                                            FT_EIS.highest_freq, 
                                            FT_EIS.lowest_freq, 
                                            FT_EIS.number_of_freqs)
                            )
            )
              
    
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
        if not Rigol:
            S, fname = create_waveform(freqs, phases, sample_freq, total_time, amax,
                        amps=amps, save=save, save_path=save_path, plot=plot)
        
        if Rigol:
            S, fname = Rigol_waveform(freqs, phases, sample_freq, total_time, amax, 
                           amps=amps)
        
        return S, fname
         


def Rigol_waveform(freqs, phases, sample_freq, total_time, amax, 
                   amps= None, save_path = r'C:\Users\BRoehrich\Desktop'):
    
    import csv   
    
    S, fname = create_waveform(freqs, phases, sample_freq, total_time,
                               amax=1, amps=amps, Rigol = True, csv = True,
                               save_path = save_path, save = False, plot=False)
    
    # RIGOL CSV Header
    headerlines = [
                ('RIGOL:DG8:CSV DATA FILE', ),
                ('TYPE:Arb', ),
                ('AMP:%s Vpp'%f'{2*amax:.4f}', ) ,
                ('PERIOD:%s S'%f'{total_time:.2E}', ),
                ('DOTS:%s'%f'{sample_freq*total_time:.0f}', ) ,
                ('MODE:INSERT', ),
                ('Sample Rate:%s'%f'{sample_freq:.2f}', ) ,
                ('AWG N:0', ),
                ('x', 'y[V]'),
                ]
    
    
    fname = 'Rigol' + fname[3:] + '.csv'
    
    try:
        if len(amps) > 1:
            fname = 'Rigol_opt' + fname[5:-4] + '.csv'
    except:
        pass
    
    file = os.path.join(save_path, fname)
    
    with open(file, 'w', newline = '\n') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        
        for line in headerlines:
            writer.writerow(line)
            
        for i in S:
            writer.writerow(('', i))
    
    print('Saved to ', os.path.join(save_path, fname))
    
    return S, fname



def update_waveform_flash_drive():
    for file in os.listdir(r'C:\Users\BRoehrich\Desktop\git\echem\HEKA\csv'):
        if file.endswith('.csv'):
            df = pd.read_csv(os.path.join(r'C:\Users\BRoehrich\Desktop\git\echem\HEKA\csv', file))
                        
            # S, fname = create_waveform(df['frequency'], df['phase'], 100000, 1, csv=True,
            #                 save = False, plot=False)
                                    
            S, fname = Rigol_waveform(df['frequency'], df['phase'], 100000, 1, amax=1, 
                           save_path=r'C:\Users\BRoehrich\Desktop\Waveforms')
            

            

