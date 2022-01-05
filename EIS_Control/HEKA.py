from process_EIS import FT_EIS
from create_waveform import waveform_from_result, Rigol_waveform
import numpy as np
import os
import time
import matplotlib.pyplot as plt
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')


filedir = r'C:\Users\BRoehrich\Desktop\HEKA python\Analysis folder'


params = {
    'out_dir': r'C:\Users\BRoehrich\Desktop\HEKA python\Analysis folder\Output',
    'csv_dir': r'C:\Users\BRoehrich\Desktop\git\echem\EIS_Control\csv',
    
    'number_of_freqs': 17,
    'highest_freq' : '17k',
    'lowest_freq' : '100',
    'filter_setting': '30k',
    'correct_filter': 1,
    'filter_correct_mode': 0,
    
    'sample_time' : 1,
    'AUTOLAB': 0,
    'autolab_current_range': 1e-6,    
    'plot_Nyquist' : 1,
    'plot_Bode' : 1,
    
    'save_Excel' : 0,
    'save_ASCII' : 0,
        }


        


def make_new_waveform(filedir, sample_freq, total_time, amax):
    
    for file in os.listdir(filedir):
        if file.endswith('.asc'):
            print('File: ', file)
            
            f = os.path.join(filedir, file)
            
            params['sample_time'] = 10
            params['correct_filter'] = True
            params['filter_correct_mode'] = False
            params['save_Excel'] = False
            params['save_ASCII'] = False
            
            FT = FT_EIS(f, params)
            
            waveform_from_result(FT, sample_freq=sample_freq, 
                                 total_time=total_time, amax=amax)
            
            
def make_new_waveform_Rigol(filedir, sample_freq, total_time, amax):
    
    for file in os.listdir(filedir):
        if file.endswith('.asc'):
            print('File: ', file)
            
            f = os.path.join(filedir, file)
            
            params['sample_time'] = 10
            params['correct_filter'] = True
            params['filter_correct_mode'] = False
            params['save_Excel'] = False
            params['save_ASCII'] = False
            
            FT = FT_EIS(f, params)
            
            S, fname = waveform_from_result(FT, sample_freq=sample_freq, 
                             total_time=total_time, amax=1, save=False,
                             Rigol = True)
            
            # Rigol_waveform(S, sample_freq, total_time, 2*amax)
            
#%%

if __name__ == '__main__':
    l = []
    n = 1
    for file in os.listdir(filedir):
        if file.endswith('.asc'):
            start_time = time.time()
            print('File %s: ' % n, file)
            
            f = os.path.join(filedir,file)
            
            FT = FT_EIS(f, params)
            
            l.append(FT)
            
            print(f'{time.time() - start_time:.2f} s\n')
            
            n += 1             
                
#%%
# Plot normalized Z vs t at select frequencies

# data = l[0].ft

# freqs = data[1]['f'].to_numpy()
# Zall = [ 
#         [data[i].iloc[j,3] for i in data] 
#         for j in range(len(freqs))
#         ]

# fig, ax = plt.subplots()

# for j in [100,1000, 4400, 17000]:
#     i = np.where(freqs == j)[0][0]
#     ax.plot(range(180)[1:], np.abs(Zall[i])[1:]/np.abs(Zall[i][1]),
#         label = f'{freqs[i]} Hz')
        
# ax.set_xlabel('Time/ s')
# ax.set_ylabel('Normalized |Z|')
# ax.legend()