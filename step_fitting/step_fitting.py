import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from step_ga import genetic_algorithm, step_func, step_func_score
import time
plt.style.use('C:/Users/BRoehrich/Desktop/git/echem/scientific.mplstyle')



file = 'C:/Users/BRoehrich/Desktop/New folder/data.asc'

# ind, time, curr, t2, v = np.loadtxt(file, unpack=True, skiprows=50, delimiter=',')

# y = curr[2000:]


arr = np.ones(1000)
for i in range(len(arr)):
    arr[i] = np.ceil(i/300)
arr += np.random.normal(loc=0, scale=0.3, size=1000)

y = arr



def initial_fit(y, n_steps):
    c, best_eval = genetic_algorithm(y, n_steps, n_pop=100,
                                        n_iter=100)
    fit = step_func(c)
    
    return c, fit


def main_fit(y, c, n_steps):
    # c: n-1 steps best candidate
    
    c[0] = (c[1][0], c[0][1])
    locs = [loc for val, loc in c]  
    
    co = []
    co.append(c)
    
    for i, (val, loc) in enumerate(c):
        if i == 0:
            continue
        
        ydata = y[c[i-1][1] : c[i][1]]
        
        best_fit, _, _ = genetic_algorithm(ydata, n_steps=1, n_pop=100,
                                       n_iter=100)
        
            
            
        
        step_loc = best_fit[1][1]
        step_val = best_fit[1][0]
        
        s_loc = c[i-1][1] + step_loc
        
        this_cand = c.copy()
        this_cand.append((step_val, s_loc))
        this_cand.sort(key=lambda x:x[1])
        
        co.append(this_cand)
    
    
    # See which new step gives the best score
    best, best_eval = co[0], step_func_score(y, co[0])
    scores = [step_func_score(y, c) for c in co]
    
        
    for i in range(len(co)):
        if scores[i] < best_eval:
            best, best_eval = co[i], scores[i]
    
        
    return best





def counter_fit(y, c):
    
    c[0] = (c[1][0], c[0][1])
    
    locs = [loc for val, loc in c]    
    
    cf = []
    cf.append(c[0])
        
    for i, (val, loc) in enumerate(c):
        if i == 0: 
            continue
        
        ydata = y[c[i-1][1] : c[i][1]]
        
        
        best_counter_fit, _, _ = genetic_algorithm(ydata, 
                                                   n_steps=1, n_pop=100,
                                                   n_iter=100)  
        
        
        
        # if best_counter_fit == 0:
        #     print(c)
        #     print(i, val, loc)
            # print(ydata)
                                                
        counter_step_loc = best_counter_fit[1][1]
        counter_step_val = best_counter_fit[1][0]
        
        cs_loc = c[i-1][1] + counter_step_loc
        cf.append((counter_step_val, cs_loc))
    
    cf.append(c[-1])
    cf.sort(key=lambda x:x[1])
        
    return cf
    
    

fits = {}
counter_fits = {}
max_steps = 20

# n_pops  = [100,200,500,1000]
# n_iters = [50,100,200]
# bestss = []
# labels=[]

# for n_iter in n_iters:
#     for n_pop in n_pops:
#         st = time.time()
#         c, best_eval, bests = genetic_algorithm(y, n_steps=10,
#                                                 n_pop=n_pop,
#                                                 n_iter=n_iter)
#         ts = time.time()-st
#         bestss.append(bests)
#         labels.append(f'N_iter={n_iter}, n_pop={n_pop}, t={ts:.2f}')
        
#         fig, ax = plt.subplots()
#         ax.plot(y)
#         ax.plot(step_func(c))
#         ax.set_title(f'N_iter={n_iter}, n_pop={n_pop}, t={ts:.2f}')
        
#         print(f'N_iter={n_iter}, n_pop={n_pop}, t={ts:.2f}')
        
# fig, ax = plt.subplots()
# for i in range(len(bestss)):
#     ax.plot(bestss[i], label=labels[i])

        


locs = []
for n_steps in range(1,5):
    if n_steps == 1:
        c, best_eval, bests = genetic_algorithm(y, n_steps=n_steps, n_pop=100,
                                            n_iter=100, force_locs=locs)
    else:
        c = main_fit(y, c, n_steps)
    
    
        
    # cf = counter_fit(y, c)
    
    fig, ax = plt.subplots()
    ax.plot(y)    
    ax.plot(step_func(c), label='fit')
    # ax.plot(step_func(cf), label='counter fit')
    ax.set_title(f'{n_steps} steps')
    
    fits[n_steps] = c
    # counter_fits[n_steps] = cf
    
    locs = [loc for (val, loc) in c if loc != 0 and loc != len(y)]
    
    # if n_steps != max_steps + 1:
        
    #     counter = counter_fit(y, fits[n_steps+1])
        
    
    
    



