import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from decay_fits import InteractivePicker

plt.ion()
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

plt.style.use('C:/Users/orozc/Google Drive (miguelorozco@ucsb.edu)/Research/Spyder/scientific.mplstyle')
data_folder = r'C:/Users/orozc/Google Drive (miguelorozco@ucsb.edu)/Research/Spyder/Run'
# plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')
# data_folder = r'C:\Users\BRoehrich\Desktop\Miguel data'


files = [os.path.join(data_folder, file) for file
         in os.listdir(data_folder) if file.endswith('.txt')]


correct_baselines = False
start_after = 100 # cut off first (n) seconds
delay = 10 # Points after "fast" spike to skip fitting on
thresh = 0.005 # Used to determine acceptable baseline "flatness"
               # Smaller = more picky, need flatter baseline to accept spike

apply_filter     = False
filter_freq      = 25

# fit = 'monoexponential'
# fit = 'monoexp-linear'
fit = 'biexponential'
# fit = 'linear'



# Define fit function and parameters    
if fit == 'biexponential':
    params=['a/ A',
            'b/ s-1',
            'c/ A',
            'd/ s-1',
            'e/ A']
    
    def exp_func(x, a,b,c,d,e):
        return a * np.exp(-b * x) + c * np.exp(-d * x) + e


elif fit == 'monoexponential':
    params = ['a/ A',
            'b/ s-1',
            'c/ A']
    
    def exp_func(x, a,b,c):
        return a * np.exp(-b * x) + c


elif fit == 'linear':
    params = ['a/ A s-1', 
              'b/ A']
    
    def exp_func(x, a, b):
        return a*x + b
    
elif fit == 'monoexp-linear':
    params = ['a/ A',
            'b/ s-1',
            'm/ A s-1',
            'c/ A']
    
    def exp_func(x, a, b, m, c):
        return a*np.exp(-b * x) + m*x + c




def store_data(ip):
        
    d = {}
    d['Index']          = []
    d['time/ s']        = []
    d.update({param: [] for param in params})
    d['Catalytic area/ C'] = []
    d['baseline slope'] = []
    d['integral/ C']    = []
    d['Reduced Chi^2']  = []
    d['File']           = []
        
            
    for i in range(len(ip.integrals)):
        d['Index'].append(i+1)
    for i in range(len(params)):
        for x in ip.save_fits[i]:
            d[params[i]].append(x)
    for x in ip.areas:
        d['Catalytic area/ C'].append(x)
    for x in ip.ms:
        d['baseline slope'].append(x)
    for x in ip.integrals:
        d['integral/ C'].append(x)
    for x in ip.chi_sqs:
        d['Reduced Chi^2'].append(x)
    for x in ip.points:
        d['time/ s'].append(x)
    
    d['File'].append(ip.file)
    d['File'].append(f'Baseline correct: {correct_baselines}')
    d['File'].append(f'Delay: {delay} pts')
    
    for key in d:
        while len(d[key]) < max([len(d[key]) for key in d]):
            d[key].append(' ')

        
    
    return d


def export(dicts, output_name):
    #dicts = list of save dicts
    dout  = {}
    for key in dicts[0].keys():
        dout[key] = list(d[key] for d in dicts)
        # Magic one liner
        # for sublist in list: for item in sublist: yield item
        dout[key] = [item for sublist in dout[key] for item in sublist]
    
    
    df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in dout.items() ]))


    writer = pd.ExcelWriter(data_folder + f'/{output_name}.xlsx', engine='xlsxwriter')
    df.to_excel(writer, index=False, header=True, startcol=0)
    writer.save()
    print(f'Saved as {data_folder}/{output_name}.xlsx')
    return df


for thresh in [0.001, 0.005, 0.01, 10]:
    ds = []
    for file in files:
        
        fig, ax = plt.subplots(figsize=(5,5), dpi=100)
        
        ip = InteractivePicker(file, fig, ax, params, exp_func, plot=False)
        
        fits = ip.analyze_file(correct_baselines, apply_filter, 
                        delay, start_after, params, exp_func, thresh)
        
        d = store_data(ip)
        ds.append(d)
        plt.close(fig)
        plt.pause(0.001)
        
    df = export(ds, f'out_{thresh}')
    
    tc1 = df['b/ s-1'].to_numpy()
    tc2 = df['d/ s-1'].to_numpy()
    
    s = [max(tc1[i], tc2[i]) for i in range(len(tc1))]    
    f = [min(tc1[i], tc2[i]) for i in range(len(tc1))]
    
    slow = []
    fast = []
    
    for i in s:
        try:
            slow.append(float(i))
        except:
            continue
    for i in f:
        try:
            fast.append(float(i))
        except:
            continue
    
    
    fig, ax = plt.subplots()
    bins = np.arange(0, 5, 0.1)
    ax.hist(fast, bins = bins, label=f'fast b, median={np.median(fast):0.2f}')
    ax.hist(slow, bins = bins, label=f'slow b, median={np.median(slow):0.2f}')
    ax.set_xlabel('b/ $s^{-1}$')
    ax.set_ylabel('Count')
    ax.legend()
    ax.set_title(f'Thresh = {thresh}, N = {len(tc1)}')
    

