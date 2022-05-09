import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


data_dir = 'C:/Users/BRoehrich/Desktop/2022-05-03/6 hr -330mV'


def createFolder(directory):
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            print ('Error: Creating directory. ' +  directory)
 
            
 
def clean_folder(folder):
    for file in os.listdir(folder):
        if file.endswith('_fit.txt'):
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
            
            
            
def extract_meisp_data(data_dir, folders):
    
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
    
    params = {'Rs': np.array([]),
              'Rct': np.array([]),
              'Cdl': np.array([]),
              'Cad': np.array([])
              }
    
    devs = {'Rs': np.array([]),
              'Rct': np.array([]),
              'Cdl': np.array([]),
              'Cad': np.array([])
              }
    
    for path in folders:
        name = path.split('/')[-1]
        
        par_file = os.path.join(path, f'{name}.par')
        dev_file = os.path.join(path, f'{name}.dev')
        
        _, _, Rs, Rct, Cdl, Cad = np.loadtxt(par_file, unpack=True)
        _, _, _, Rs_dev, Rct_dev, Cdl_dev, Cad_dev = np.loadtxt(dev_file, 
                                                          unpack=True)
                                                         
        params['Rs'] = np.concatenate((params['Rs'], Rs))
        params['Rct'] = np.concatenate((params['Rct'], Rct))
        params['Cdl'] = np.concatenate((params['Cdl'], Cdl))
        params['Cad'] = np.concatenate((params['Cad'], Cad))
        devs['Rs'] = np.concatenate((devs['Rs'], Rs_dev))
        devs['Rct'] = np.concatenate((devs['Rct'], Rct_dev))
        devs['Cdl'] = np.concatenate((devs['Cdl'], Cdl_dev))
        devs['Cad'] = np.concatenate((devs['Cad'], Cad_dev))
    
    params['ket'] = 1/(2*params['Rct']*params['Cad'])    
    
    return times, params, devs



def save_as_excel(times, params, devs):
    df = pd.DataFrame([times, 
                   *[val for key, val in params.items()], 
                   *[val for key, val in devs.items()]])
    df = df.T
    df = df.rename({0: 'time/s',
           1: 'Rs',
           2: 'Rct',
           3: 'Cdl',
           4: 'Cad',
           5: 'ket',
           6: 'Rs_dev',
           7: 'Rct_dev',
           8: 'Cdl_dev',
           9: 'Cad_dev'}, axis=1)
    
    writer = pd.ExcelWriter(data_dir + '/meisp_fits.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {data_dir}/meisp_fits.xlsx')
    
    return f'{data_dir}/meisp_fits.xlsx'
    
    


# folder_numbers = make_meisp_folders(data_dir)
# print(folder_numbers)
folder_numbers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

times, params, devs = extract_meisp_data(data_dir, folder_numbers)

output = save_as_excel(times, params, devs)







#%%


# # Plot params vs t & normalized params vs t for invivo measurement
file = output

df = pd.read_excel(file)


from scipy import signal
fig, ax = plt.subplots()
y = signal.savgol_filter(np.array(df['ket']), 21, 1, mode='nearest')
ax.plot(df['time/s']/60, df['ket'], color = colors[0])
# ax.plot(df['time/s']/60, y, color = colors[1])
ax.set_xlabel('Time/ min')
ax.set_xticks([0,30,60,90,120])
ax.set_ylabel('$k_{et}$/ $s^{-1}$')

fig, ax = plt.subplots()
for param in ['Rs', 'Rct', 'Cdl', 'Cad']:
    ax.plot(df['time/s']/60, df[param]/df[param][0], label=param)

ax.set_xlabel('Time/ min')
ax.set_ylabel('Normalized parameter')
ax.set_xticks([0,30,60,90, 120])
# ax.set_xlim(ax.get_xlim()[0], 250)
ymin, ymax = ax.get_ylim()
ax.set_ylim(0.3, ymax)
ax.legend()



            
