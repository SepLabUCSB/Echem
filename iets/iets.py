import pyvisa
import time
import os

from funcs import init_recording, record_current


scope   = 'USB0::0xF4ED::0xEE3A::SDS1EDED5R0471::INSTR'
autolab_i_range    = 10e-6 # A/V
lockin_sensitivity = 500e-3 # V


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
    params = init_recording(inst)
    
    # if os.path.exists(trig_file):
    #     os.remove(trig_file)
        
    j = 0
    
    while j < n:
        
        while os.path.exists(trig_file) == False:
            time.sleep(0.001)
            
        # record single (v, i)
        i = record_current(inst, 0, params, autolab_i_range,
                       lockin_sensitivity)
        
        j_, v = get_v(trig_file)
        
        l.append((v,i))
        if SAVE:
            append_save(save_file, (v,i))
        print(j, v, i)
        j += 1
        # os.remove(trig_file)
    
    
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
            
# l = record_sweep(scope, 51, save_path)
d = record_multiple_sweeps(scope, 3, 51, save_path)
            
