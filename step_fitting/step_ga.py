import numpy as np
import matplotlib.pyplot as plt
from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Array
plt.style.use('Z:/Projects/Brian/scientific.mplstyle')
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

'''
https://machinelearningmastery.com/simple-genetic-algorithm-from-scratch-in-python/
'''

def least_sq(data, fit):
    return np.sum((data-fit)**2)/np.mean(data)**2



def step_func(c):
    # c = [(0, 100),
    #      (1, 400),
    #      (2,1000)]    
        
    
    arr = np.array([])
        
    # if not any(pt[1] == 0 for pt in c):
    #     c.append((c[0][0], 0))
    c.sort(key=lambda x:x[1])
        
    for i, (val, end) in enumerate(c):
        if i > 0:
            last_idx = c[i-1][1]
            new_arr = val*np.ones(end - last_idx)
            arr = np.hstack((arr, new_arr))
        
    return arr


def step_func_score(y, candidate):
    arr = step_func(candidate)
    return least_sq(y, arr)


def random(lower, upper):
    return (upper - lower) * np.random.random_sample() + lower


def selection(pop, scores, k=3):
    '''
    Parameters
    ----------
    pop : list of candidates.
    scores : list of candidate scores.
    k : int.
    
    Choose k random individuals, choose parent as fittest
    of those individuals

    '''
    # Choose k random individuals from population
    # Select parent as fittest of those individuals
    selection_ix = np.random.randint(len(pop))
    
    for ix in np.random.randint(0, len(pop), k-1):
        if scores[ix] < scores[selection_ix]:
            selection_ix = ix
    
    return pop[selection_ix]


def crossover(p1, p2, r_cross, force_locs):
    '''
    Create children c1 and c2 which are linear combinations of
    parents p1 and p2.
    '''
    #children are copies of parents by default
    c1, c2 = p1.copy(), p2.copy()
    
    if np.random.rand() < r_cross:
        
        frac = np.random.rand()
        
                
        for i in range(len(c1)):
            
            
            loc1 = int(frac*p1[i][1] + (1-frac)*p2[i][1])
            loc2 = int(frac*p2[i][1] + (1-frac)*p1[i][1])
            
            if p1[i][1] in force_locs:
                loc1 = p1[i][1]
            if p2[i][1] in force_locs:
                loc2 = p2[i][1]
            
            
            
            c1[i] = (frac*p1[i][0] + (1-frac)*p2[i][0], 
                     loc1)
            
            c2[i] = (frac*p2[i][0] + (1-frac)*p1[i][0], 
                     loc2)
            
    
    c1.sort(key=lambda x:x[1])
    c2.sort(key=lambda x:x[1])
            
        
    return [c1, c2]




def mutation(candidate, bounds, length, r_mut, force_locs):
    '''
    Mutate candidate with r_mut chance
    
    Mutation performed on individual parameters,
    
    randomizes parameter value within set bounds
    '''  
    
        
    for i, (val, loc) in enumerate(candidate):
        
        locs = [loc for (val, loc) in candidate]
        
        if (np.random.rand() < r_mut and loc != 0):
            val = random(bounds[0], bounds[1])
        
        if (np.random.rand() < r_mut and loc != 0 and loc != length
            and loc not in force_locs):
            loc = np.random.randint(1, length)
            
            bad_locs = [loc for loc in (loc-1, loc, loc+1)
                        if loc in locs]
            
            if len(bad_locs) != 0:
                
                if loc > length:
                    loc = np.random.randint(1, length)
                # while any((loc-1, loc, loc+1)) in locs:
                #     print(loc, locs)
                #     loc += 1
                # print('changed loc: ', loc, locs)
        
        candidate[i] = (val, loc)
    
    if len(np.unique([loc for (val, loc) in candidate])) != len(candidate):
        candidate = mutation(candidate, bounds, length, 1, [])
    
    # candidate.sort(key=lambda x:x[1])
    # candidate[0] = (candidate[1][0], candidate[0][1])
    
    return candidate
    

def worker(l, i, candidate, y, objective):
    score = objective(y, candidate)
    l[i] = score
    



def genetic_algorithm(y, n_steps, 
                       objective=step_func_score, n_pop=100, 
                      n_iter=100, r_cross=0.9, r_mut=0.25, ax=None,
                      starting_guess = None, force_locs = []):
    
    
    #generate initial pop
    pop = []
    bests = []
    
    
    bounds = [min(y), max(y)]
    length = len(y)
    
    for i in range(n_pop):
 
        c = []
        locs = []
        
        if len(force_locs) != 0:
            for loc in force_locs:
                locs.append(loc)
                c.append(
                    (random(bounds[0], bounds[1]), loc)
                    )
        
        for _ in range(n_steps - len(c)):
            
            loc = np.random.randint(1, length)
            
            if loc in locs:
                while (loc or loc-1 or loc+1) in locs:
                    loc = np.random.randint(1, length)
            
            c.append(
                (random(bounds[0], bounds[1]), loc)
                )
        
                
        # Endpoint
        c.append(
            (random(bounds[0], bounds[1]), 
             length)
                )
        
        c.sort(key=lambda x:x[1])
        
        # Starting point
        c.append((c[0][0], 0))
        c.sort(key=lambda x:x[1])
                
        pop.append(c)
                
    
    
    if starting_guess:
        for i in range(n_pop):
            if i < int(0.1*n_pop):
                # Make 10% of the population the starting guess
                pop[i] = starting_guess
   

      
        
    # initialize best candidate
    best, best_eval = pop[0], objective(y, pop[0])
    
    # if force_val:
    #     for i in range(n_pop):
    #         print(i, pop[i])
    #         print('')
    
    # iterate through generations
    for gen in range(n_iter):
        #check for new best solution
        scores = [objective(y, c) for c in pop]
        
        for i in range(n_pop):
            if scores[i] < best_eval:
                best, best_eval = pop[i], scores[i]
                
                # print("Gen %d, new best = %.3f" % (gen, scores[i]))

        if (ax and gen%10 == 0):
            # plot iteration every 10 generations
            colors = plt.cm.magma(np.linspace(0,0.8,n_iter))
            
            fit_arr = step_func(best)
            
            ax.plot(fit_arr,color = colors[gen], 
                    alpha = (0.5*gen/n_iter)+0.2)

        
        #select parents
        selected = [selection(pop, scores) for _ in range(n_pop)]
        
        #create next generation
        children = list()
        for i in range(0, n_pop, 2):
            #get pairs of selected parents
            p1, p2 = selected[i], selected[i+1]
            
            #crossover and mutate
            for c in crossover(p1, p2, r_cross, force_locs):
                c = mutation(c, bounds, length, r_mut, force_locs)
                children.append(c)
        
        # for i in range(len(children)):
        #     # Duplicate 10% of population as previous best candidate
        #     if i < int(0.1*n_pop):
        #         children[i] = best

        bests.append(best_eval)
        
        #replace population
        pop = children
        
    
    # fig, ax2 = plt.subplots()
    # ax2.plot(bests)
    
    return [best, best_eval, bests]



if __name__ == '__main__':
    file = 'C:/Users/BRoehrich/Desktop/New folder/data.asc'
    
    ind, time, curr, t2, v = np.loadtxt(file, unpack=True, skiprows=50, delimiter=',')
    
    
    x = np.linspace(0,1, 1000)
    y1 = np.zeros(300)
    y2 = np.ones(300)
    y3 = 1.5*np.ones(400)
    y = np.concatenate((y1, y2, y3)) + np.random.rand(1000) - 0.5
    
    # y = curr[2000:]
    
    fig, ax = plt.subplots()
    ax.plot(y)
    
    best, best_eval, bests = genetic_algorithm(y, n_steps=2, ax=ax)


