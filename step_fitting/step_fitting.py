import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from step_ga import genetic_algorithm, step_func, step_func_score
import time
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')

## Adapted from Nature 442, 709–712 (2006). https://doi.org/10.1038/nature04928 ##

# file = 'C:/Users/BRoehrich/Desktop/New folder/data.asc'
file = 'C:/Users/BRoehrich/Desktop/Sofia data/3 CVA 2 mM FeMeOH half mM KCl -400 mV 20 min 50 Hz CUME26 flow rate 20uLhr to 020uLhr_C02.txt'
time, curr = np.loadtxt(file, unpack=True, skiprows=1, delimiter='\t')
y = curr[5000:]
from scipy.signal import savgol_filter
filt = savgol_filter(y, 101, 1)
filt2 = np.average(filt.reshape((-1, 5)), axis=1)
plt.plot(filt2)
y = filt2
# ind, time, curr, t2, v = np.loadtxt(file, unpack=True, skiprows=50, delimiter=',')

# y = curr[2000:]
#%%


# arr = np.ones(1000)
# for i in range(len(arr)):
#     arr[i] = np.ceil(i/100)
# arr += np.random.normal(loc=0, scale=0.6, size=1000)

# y = arr



def initial_fit(y, n_steps):
    c, best_eval = genetic_algorithm(y, n_steps, n_pop=100,
                                        n_iter=100)
    fit = step_func(c)
    
    return c, fit


def main_fit(y, c, n_steps):
    # c: n-1 steps best candidate

    co = []
    
    for i, (loc, val) in enumerate(c):
        if i == 0:
            continue
        
        ydata = y[c[i-1][0] : c[i][0]]
        
        if len(ydata) < 5:
            continue
        
        best_fit, _, _ = genetic_algorithm(ydata, n_steps=1, n_pop=100,
                                       n_iter=25)
        
        
        step_loc = best_fit[1][0]
        step_val = best_fit[2][1]
        
        s_loc = c[i-1][0] + step_loc
        
        this_cand = c.copy()
        this_cand.append((s_loc, step_val))
        this_cand.sort()
        
        co.append(this_cand)
    
    
    # Reevaluate median values
    for c in co:
        for i in range(1, len(c)):
            med = np.median(y[ c[i-1][0]:c[i][0] ])
            c[i] = (c[i][0], med)
        c[0] = (c[0][0], c[1][1])
    
    
    # See which new step gives the best score
    best, best_eval = co[0], step_func_score(y, co[0])
    scores = [step_func_score(y, c) for c in co]
        
    for i in range(len(co)):
        if scores[i] < best_eval:
            best, best_eval = co[i], scores[i]
    
        
    return best, best_eval





def counter_fit(y, c):
    
    cf = []
    cf.append(c[0])
        
    for i, (loc, val) in enumerate(c):
        if i == 0: 
            continue
        
        ydata = y[c[i-1][0] : c[i][0]]
        
        if len(ydata) < 5:
            continue
        
        
        best_counter_fit, _, _ = genetic_algorithm(ydata, 
                                                   n_steps=1, n_pop=100,
                                                   n_iter=25)  
        
        
        
        # if best_counter_fit == 0:
        #     print(c)
        #     print(i, val, loc)
            # print(ydata)
                                                
        counter_step_loc = best_counter_fit[1][0]
        counter_step_val = best_counter_fit[2][1]
        
        cs_loc = c[i-1][0] + counter_step_loc
           
        cf.append((cs_loc, counter_step_val))
        
        
        
    cf.append(c[-1])
    # for i in range(1, len(c)):
    #     cf.pop(i)
    cf.sort()
    
    # Pick counter fit candidate with the best n steps
    best      = None
    best_eval = 1e20
    # print('')
    # print('c: ',  [loc for (loc, val) in c])
    # print('cf: ', [loc for (loc, val) in cf])
    for i in range(1, len(cf)-1):
        cf_cand = cf.copy()
        
        cf_cand.pop(i)
        
        for i in range(1, len(cf_cand)):
            med = np.median(y[ cf_cand[i-1][0]:cf_cand[i][0] ])
            cf_cand[i] = (cf_cand[i][0], med)
        
        cf_cand[0] = (cf_cand[0][0], cf_cand[1][1])
        score = step_func_score(y, cf_cand)
        
        # print(score, [loc for (loc, val) in cf_cand])
        
        if score < best_eval:
            best      = cf_cand
            best_eval = score
        
    return best, best_eval
    


def S(fit_score, cf_score):
    return cf_score/fit_score
    


def refine_locs(y, c, width):
    
    # Look within +- width/2 points for the largest step
    
    locs = [loc for i, (loc, val) in enumerate(c)]
    tol = int(width/10)
    
    for i in range(1, len(c)):
        loc = c[i][0]
        if any([x in locs for x in range(loc-tol, loc+tol+1) if x != loc]):
            print(f'Removing point at {loc}')
            c.pop[i]
    
    for i, (loc, val) in enumerate(c):
        
        if i == 0 or i == len(c)-1:
            continue
        
        ydata = y[int(loc-width/2) : int(loc+width/2)]
         
        best_fit, _, _ = genetic_algorithm(ydata, n_steps=1, n_pop=100,
                                           n_iter=25)
        
        
        new_loc = int(loc-width/2) + best_fit[1][0]
        
        print(loc, new_loc)
        
        c[i] = (new_loc, val)
    
    
    # Recalculate median values
    
    for i in range(1, len(c)):
         if i > 0:
             med = med = np.median(y[ c[i-1][0]:c[i][0] ])
             c[i] = (c[i][0], med)
    
    c[0] = (c[0][0], c[1][1])
    
    return c



fits = {}
counter_fits = {}
Ss    = []
steps = []
max_steps = 20


        
for n_steps in range(1,20):
    if n_steps == 1:
        c, c_score, bests = genetic_algorithm(y, n_steps=n_steps, n_pop=100,
                                            n_iter=25)
    else:
        c, c_score = main_fit(y, c, n_steps)
    
    
        
    cf, cf_score = counter_fit(y, c)
    
    fig, ax = plt.subplots()
    ax.plot(y)    
    ax.plot(step_func(c), label='fit')
    ax.plot(step_func(cf), label='counter fit')
    ax.set_title(f'{n_steps} steps')
    
    fits[n_steps] = (c, c_score)
    counter_fits[n_steps] = (cf, cf_score)
    
    locs = [loc for (val, loc) in c if loc != 0 and loc != len(y)]
    
    Ss.append(S(c_score, cf_score))
    steps.append(n_steps)
    
    print(f'{n_steps} steps: S = {S(c_score, cf_score)}')


fig, ax = plt.subplots()
ax.plot(steps, Ss)
ax.set_xlabel('# of steps')
ax.set_ylabel('Score')    

print('')
print(f'I think there are {steps[np.where(Ss == max(Ss))[0][0]]} steps')

idx = np.where(Ss == max(Ss[1:]))[0][0]+1
refined_c = refine_locs(y, fits[idx][0], len(y)/10)

fig, ax = plt.subplots()
ax.plot(y)
ax.plot(step_func(fits[idx][0]))
ax.plot(step_func(refined_c))
        
    
    
    



