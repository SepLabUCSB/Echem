import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


data_dir = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-04-27\6hr -330mV'


def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
 
            
 
def clean_folder(folder):
    for file in os.listdir(folder):
        if not file.endswith('s.txt'):
            os.remove(os.path.join(folder,file))
            
            

def make_meisp_folders(data_dir):
    i = 0            
    for file in os.listdir(data_dir):
        if len(file) == 9 and file.endswith('s.txt'):
            i += 1
    
    no_of_folders = np.ceil(i/450).astype(int)
    
    for i in range(no_of_folders):
        new_folder = os.path.join(data_dir, str(i))
        createFolder(new_folder)
        print(new_folder)
        idx = 450*i
        idy = 450*(i+1)
        for j in range(idx, idy):
            try:
                src = os.path.join(data_dir, f'{j:04}s.txt')
                dst = os.path.join(new_folder, f'{j:04}s.txt')
                shutil.copyfile(src, dst)
            except FileNotFoundError:
                continue
    
    return [i for i in range(no_of_folders)]
            
            
            
def extract_meisp_data(data_dir, folders, circuit_elements=None):
    
    time_file = data_dir + '/0000_time_list.txt'
    times = []
    with open(time_file, 'r') as f:
        for line in f:
            try:
                time = float(line)
                times.append(time)
            except:
                pass
    
    folders = [data_dir + f'/{str(folder)}' for folder in folders]
    
    
    if not circuit_elements:
        circuit_elements = ['Rs', 'Rct', 'Cdl', 'Cad']
    
    params = {elem: np.array([]) for elem in circuit_elements}
    devs   = {elem: np.array([]) for elem in circuit_elements}
    
    
    
    for path in folders:
        name = path.split('/')[-1]
        
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



def save_as_excel(times, params, devs):
    df = pd.DataFrame([times, 
                   *[val for key, val in params.items()], 
                   *[val for key, val in devs.items()]])
    df = df.T

    names = ['time/s'] + list(params) + [f'{dev}_dev' for dev in devs]
    df = df.rename({i: names[i] for i in range(len(names))}, axis=1)
            
    
    writer = pd.ExcelWriter(data_dir + '/meisp_fits.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {data_dir}/meisp_fits.xlsx')
    
    return f'{data_dir}/meisp_fits.xlsx'
    
    


# folder_numbers = make_meisp_folders(data_dir)
# print(folder_numbers)
folder_numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

# Pass list of circuit elements in the order they appear in MEISP
circuit_elements = ['Rs', 'Rct', 'Cad', 'phi', 'Cdl']
times, params, devs = extract_meisp_data(data_dir, folder_numbers,
                                          circuit_elements = circuit_elements)

output = save_as_excel(times, params, devs)







#%%


# # Plot params vs t & normalized params vs t for invivo measurement
file = output

df = pd.read_excel(file)


# df = df[df['ket'] > 70]
# df = df[df['Cad'] > 0]


from scipy import signal
fig, ax = plt.subplots()
y = signal.savgol_filter(np.array(df['ket']), 21, 1, mode='nearest')
ax.plot(df['time/s']/60, df['ket'], color = colors[0])
# ax.plot(df['time/s']/60, y, color = colors[1])
ax.set_xlabel('Time/ min')
ax.set_xticks([0,30,60,90])
ax.set_ylabel('$k_{et}$/ $s^{-1}$')

fig, ax = plt.subplots()
for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
    ax.plot(df['time/s']/60, df[param]/df[param][0], label=param)

ax.set_xlabel('Time/ min')
ax.set_ylabel('Normalized parameter')
ax.set_xticks([0,30,60,90])
# ax.set_xlim(ax.get_xlim()[0], 250)
ymin, ymax = ax.get_ylim()
ax.set_ylim(0.3, ymax)
ax.legend()



            
