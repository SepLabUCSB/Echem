import pyvisa
import time
import os
import matplotlib.pyplot as plt
import numpy as np

from funcs import init_recording, record_current, record_current_srs, record_current_sr

plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')


scope   = 'USB0::0xF4ED::0xEE3A::SDS1EDED5R0471::INSTR'
srs     = 'GPIB0::8::INSTR'
sr      = 'GPIB0::16::INSTR'
autolab_i_range    = 1e-3 # A
lockin_sensitivity = 10e-3 # V


trig_file = r'C:\Users\BRoehrich\Desktop\git\echem\iets\update.txt'
save_path = r'C:\Users\BRoehrich\Desktop\iets'


def record_multiple_sweeps(scope, n_sweeps, n_steps, save_path):
    
    # Create folder to hold sweeps
    folder = input('Folder name >> ')
    folder = os.path.join(save_path, folder)
    createFolder(folder)
    
    d = {}
    
    for n in range(n_sweeps):
        l = record_sweep(scope, n_steps, folder, SAVE=True,
                         save_file=f'{n}')
        d[n] = l

    return d


def record_sweep(scope, n, save_path, ask_inputs=False, SAVE=False,
                 save_file=''):
    
    
    if ask_inputs:
        SAVE = input('Save? (0/1) >> ')
        
        if SAVE == '1': SAVE = True
        else: SAVE = False
        
    if SAVE:
        if ask_inputs:
            save_file = input('Save name? >> ')
        save_file = os.path.join(save_path, f'{save_file}.txt')
        
        with open(save_file, 'w') as f:
            f.write('DC V, AC rms A\n')
            f.close()
    
    
    l = []
    
    # connect to scope
    rm = pyvisa.ResourceManager()
    inst = rm.open_resource(scope)
    
    
    # initialize scope params
    if scope.startswith('USB'):
        params = init_recording(inst)
    # elif scope.startswith('GPIB'):
    #     inst.write('AGAN')
    
    if scope == 'GPIB0::16::INSTR':
        sens = inst.query('SENS')
        sens = int(sens)
    
    if os.path.exists(trig_file):
        os.remove(trig_file)
        
    j = 0
    
    while j < n:
        
        # if (j % 50 == 0) and scope.startswith('GPIB'):
        #     print('Previous sensitivity ', inst.query('SENS?'))
        #     inst.write('AGAN')
            
        
        while os.path.exists(trig_file) == False:
            time.sleep(0.001)
            
        # record single (v, i)
        if scope.startswith('USB'):
            i = record_current(inst, 0, params, autolab_i_range,
                           lockin_sensitivity)
            
        elif scope == 'GPIB0::8::INSTR':
            i = record_current_srs(inst, autolab_i_range, 0.3)
            
        elif scope == 'GPIB0::16::INSTR':
            i = record_current_sr(inst, sens, autolab_i_range, 0.4)
        
        j, v = get_v(trig_file)
        
        l.append((v,i))
        if SAVE:
            append_save(save_file, (v,i))
        print(j, v, i)
        
        os.remove(trig_file)
    
    
    return l


def get_v(file):
    with open(file, 'r') as f:
        i = 0
        for line in f:
            if i == 2:
                num, v = line.strip('\n').split('\t')
                return int(num), float(v)
            i += 1
            
            
def append_save(save_file, dtup):
    v, i = dtup
    with open(save_file, 'a') as f:
        f.write(f'{v},{i}\n')
        f.close()


def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)            

            
# l = record_sweep(sr, 201, save_path, ask_inputs=True, SAVE=True)
d = record_multiple_sweeps(sr, 100, 201, save_path)




# fig, ax = plt.subplots()
# for n in d:
#     l = d[n]
#     ii = []
#     vv = []
#     for (v, i) in l:
#         ii.append(i)
#         vv.append(v)
    
#     ax.plot(vv, ii)

# v = np.average([[v for (v,i) in d[n]] for n in d], axis=0)
# i = np.average([[i for (v,i) in d[n]] for n in d], axis=0)
            
