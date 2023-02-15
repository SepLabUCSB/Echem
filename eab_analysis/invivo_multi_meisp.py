import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
import shutil
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

mpl_is_inline = 'inline' in matplotlib.get_backend()
if not mpl_is_inline:
    plt.rcParams['figure.figsize'] = (5,5)
    plt.rcParams['figure.dpi'] = 100


data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-09-01\-340mV vanco'


def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
 
            
 
def clean_folder(folder):
    for file in os.listdir(folder):
        if not file.endswith('.txt'):
            os.remove(os.path.join(folder,file))
        if file.endswith('fit.txt'):
            os.remove(os.path.join(folder, file))
            
            

def make_meisp_folders(data_dir, probe=''):        
    
    files   = []
    folders = []
    
    for file in os.listdir(data_dir) :
        if file.startswith(probe) and file.endswith('.txt'):
            num = int(file.strip('.txt').split('_')[-1])
            files.append((file, num))
    
    files.sort(key=lambda x:x[1])
    
    files, nums = zip(*files)
    
    no_of_folders = np.ceil(len(files)/400).astype(int)
    
    for i in range(no_of_folders):
        new_folder = os.path.join(data_dir, probe+str(i))
        createFolder(new_folder)
        print(new_folder)
        folders.append(new_folder)
        idx = 400*i
        idy = 400*(i+1)
        for j in nums[idx:idy]:
            try:
                src = os.path.join(data_dir, f'{probe}_{j}.txt')
                dst = os.path.join(new_folder, f'{probe}_{j:05}.txt')
                if not os.path.exists(dst):
                    shutil.copyfile(src, dst)
            except FileNotFoundError:
                print(f'File not found: {src}')
                continue
    
    return folders, nums
            
            
            
def extract_meisp_data(data_dir, folders, time_nums,
                       probe='', circuit_elements=None):
    
    time_file = data_dir + '/0000_time_list.txt'
    times = np.loadtxt(time_file)   
    times = times[tuple([time_nums])] # np weirdness, get times from indices
    
    if not circuit_elements:
        circuit_elements = ['Rs', 'Rct', 'Cdl', 'Cad']
    
    params = {elem: np.array([]) for elem in circuit_elements}
    devs   = {elem: np.array([]) for elem in circuit_elements}
    
    
    
    for path in folders:
        name = path.split('\\')[-1]
        
        par_file = os.path.join(path, f'{name}.par')
        dev_file = os.path.join(path, f'{name}.dev')
        
        
        ps = pd.read_fwf(par_file, names=['ind', '_', *circuit_elements])
        ds = pd.read_fwf(dev_file, names=['ind', '_', 's', *circuit_elements])
               
        
        for param in circuit_elements:
            params[param] = np.concatenate((params[param], 
                                           np.array(ps[param])))
            devs[param] = np.concatenate((devs[param], 
                                         np.array(ds[param])))
            
             
    params['ket'] = 1/(2*params['Rct']*params['Cad'])    
    
    return times, params, devs




    


def save_as_excel(times, params, devs, probe):
    df = pd.DataFrame([times, 
                   *[val for key, val in params.items()], 
                   *[val for key, val in devs.items()]])
    df = df.T

    names = ['time/s'] + list(params) + [f'{dev}_dev' for dev in devs]
    df = df.rename({i: names[i] for i in range(len(names))}, axis=1)
    
    df = df[df['Cdl'] < 2e-6]
            
    
    writer = pd.ExcelWriter(data_dir + f'/{probe}_meisp_fits.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {data_dir}/{probe}_meisp_fits.xlsx')
    
    return f'{data_dir}/{probe}_meisp_fits.xlsx'



def plot_output(file, name):   
    df = pd.read_excel(file)
        
    
    from scipy import signal
    fig, ax = plt.subplots(dpi=100)
    y = signal.savgol_filter(np.array(df['ket']), 51, 1, mode='nearest')
    ax.plot(df['time/s']/60, df['ket'], color = colors[0])
    ax.plot(df['time/s']/60, y, color = colors[1])
    ax.set_xlabel('Time/ min')
    ax.set_xticks(np.arange(0, ax.get_xlim()[1], 30))
    ax.set_ylim(0, 50)
    ax.set_ylabel('$k_{et}$/ $s^{-1}$')
    ax.set_title(name)
    fig.tight_layout()
    
    fig, ax = plt.subplots()
    for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
        ax.plot(df['time/s']/60, df[param]/df[param][0], label=param)
    
    ax.set_xlabel('Time/ min')
    ax.set_ylabel('Normalized parameter')
    ax.set_xticks(np.arange(0, ax.get_xlim()[1], 30))
    # ax.set_xlim(ax.get_xlim()[0], 250)
    ymin, ymax = ax.get_ylim()
    # ax.set_ylim(0.3, ymax)
    ax.set_title(name)
    ax.legend()
    fig.tight_layout()

#%%


rvein_folders, rvein_nums = make_meisp_folders(data_dir, probe='rvein')
lvein_folders, lvein_nums = make_meisp_folders(data_dir, probe='lvein')


#%%

# Pass list of circuit elements in the order they appear in MEISP
circuit_elements = ['Rs', 'Rct', 'Cad', 'phi', 'Cdl']


times, params, devs = extract_meisp_data(data_dir, rvein_folders, rvein_nums, 'Right vein',
                                          circuit_elements = circuit_elements)
rvein_output = save_as_excel(times, params, devs, 'rvein')




times, params, devs = extract_meisp_data(data_dir, lvein_folders, lvein_nums, 'Left vein',
                                          circuit_elements = circuit_elements)
lvein_output = save_as_excel(times, params, devs, 'lvein')



plot_output(rvein_output, 'Right Vein')
plot_output(lvein_output, 'Left Vein')





            
