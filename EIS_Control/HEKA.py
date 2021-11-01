from process_EIS import FT_EIS
from create_waveform import waveform_from_result, Rigol_waveform
import os
import matplotlib.pyplot as plt
plt.style.use('Z:\Projects\Brian\scientific.mplstyle')


filedir = r'C:\Users\BRoehrich\Desktop\HEKA python\Analysis folder'


params = {
    'out_dir': r'C:\Users\BRoehrich\Desktop\HEKA python\Analysis folder',
    'csv_dir': r'C:\Users\BRoehrich\Desktop\git\echem\HEKA\csv',
    
    'number_of_freqs': 31,
    'highest_freq' : '1k',
    'lowest_freq' : '1',
    'filter_setting': '30k',
    'correct_filter': True,
    'filter_correct_mode': False,
    
    'sample_time' : 10,
    'AUTOLAB': 1,
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

# if __name__ == '__main__':
#     l = []
#     n = 1
#     for file in os.listdir(filedir):
#         if file.endswith('.asc'):
#             print('File %s: ' % n, file)
            
#             f = os.path.join(filedir,file)
            
#             FT = FT_EIS(f, params)
            
#             l.append(FT)
            
#             n += 1             
                

