import numpy as np
import matplotlib.pyplot as plt
from EIS_Fit import circuits
import pandas as pd
# plt.style.use('Z:/Projects/Brian/scientific.mplstyle')
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

'''
https://machinelearningmastery.com/simple-genetic-algorithm-from-scratch-in-python/
'''


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


def crossover(p1, p2, r_cross):
    '''
    Create children c1 and c2 which are linear combinations of
    parents p1 and p2.
    '''
    #children are copies of parents by default
    c1, c2 = p1.copy(), p2.copy()
    
    if np.random.rand() < r_cross:
        
        frac = np.random.rand()
        
        for p in c1:
            c1[p] = frac*p1[p] + (1-frac)*p2[p]
            c2[p] = frac*p2[p] + (1-frac)*p1[p]
            
        
    return [c1, c2]




def mutation(candidate, bounds, r_mut):
    '''
    Mutate candidate with r_mut chance
    
    Mutation performed on individual parameters,
    
    randomizes parameter value within set bounds
    '''    
    for param, val in candidate.items():
        
        if np.random.rand() < r_mut:
            # Randomize parameter within set bounds
            if 'n' in param:
                # Don't pick CPE phase by log value
                if bounds[param][0] == bounds[param][1]:
                    # Force CPE phase to a certain value
                    candidate[param] = bounds[param][0]
                else:
                    candidate[param] = np.random.uniform(bounds[param][0],
                                                         bounds[param][1])

                
            else:
                # Uniform spacing on log scale
                logval = np.random.uniform(np.log10(bounds[param][0]), 
                                            np.log10(bounds[param][1]))

                candidate[param] = 10**logval

    return candidate
    




def genetic_algorithm(freqs, Z, bounds, circuit, 
                      objective=circuits.leastsq_errorfunc, n_pop=500, 
                      n_iter=200, r_cross=0.9, r_mut=0.25, ax=None,
                      starting_guess = None):
    '''
    Main function for performing genetic algorithm optimization.

    Parameters
    ----------
    freqs : array of frequencies
    
    Z : array, in form re + j*im
    
    bounds : dict, in form {'R':[low, high],}
    
    circuit : string corresponding to preset circuit in circuits.py
    
    objective : func, optional
        Objective function to minimize. 
        The default is circuits.leastsq_errorfunc.
 
    n_pop : int, optional
        Size of the population.
        
    n_iter : int, optional
        Number of generations to iterate over.
        
        
    r_cross : float between 0 and 1, optional
        Crossover rate. Children are identical to their parents
        unless rand() < r_cross.
        
    r_mut : float between 0 and 1, optional
        Mutation rate.
    
    ax: matplotlib.Axes
        axis to iteratively plot best performers onto
        
    starting_guess: dict
        dict of parameters to include as an
        individual in the starting population (adds
        1% of the starting pop as starting_guess)
    
    Returns
    -------
    list
        dict of best-fit parameters,
        score of best individual.

    '''

    #generate initial pop
    pop = []

    for i in range(n_pop):
        c = {}
        for param, lims in bounds.items():
            if 'n' in param:
                # Don't pick CPE phase by log value
                if lims[0] == lims[1]:
                    # Force CPE phase to a certain value
                    c[param] = lims[0]
                else:
                    c[param] = np.random.uniform(lims[0], lims[1])

                
            else:
                # Uniform spacing on log scale
                logval = np.random.uniform(np.log10(lims[0]), 
                                            np.log10(lims[1]))

                c[param] = 10**logval
                # c[param] = np.random.uniform(lims[0],lims[1])

        # pop.append(c)
        
        if type(c) == dict:
            pop.append(c)
    

    
    if starting_guess:
        for i in range(n_pop):
            if i < int(0.1*n_pop):
                # Make 10% of the population the starting guess
                pop[i] = starting_guess
                
        
    # initialize best candidate
    best, best_eval = pop[0], objective(freqs, Z, pop[0], circuit)

    
    # iterate through generations
    for gen in range(n_iter):
        
        #check for new best solution
        scores = [objective(freqs, Z, c, circuit) for c in pop] 
        
        for i in range(n_pop):
            if scores[i] < best_eval:
                best, best_eval = pop[i], scores[i]
                # print("Gen %d, new best = %.3f" % (gen, scores[i]))

        if ax:
            # plot iteration every 10 generations
            if gen%10 == 0:
                colors = plt.cm.magma(np.linspace(0,0.8,n_iter))
                fit_Z = circuits.chosen_circuit(circuit, freqs, best)
                ax.plot(np.real(fit_Z)/1e6, -np.imag(fit_Z)/1e6, 
                        color = colors[gen], alpha = (0.5*gen/n_iter)+0.2)
                ax.set_xlabel("Z'/ M$\Omega$")
                ax.set_ylabel("Z''/ M$\Omega$")
                ax.set_ylim(ax.get_xlim())

        
        #select parents
        selected = [selection(pop, scores) for _ in range(n_pop)]
        
        #create next generation
        children = list()
        for i in range(0, n_pop, 2):
            #get pairs of selected parents
            p1, p2 = selected[i], selected[i+1]
            
            #crossover and mutate
            for c in crossover(p1, p2, r_cross):
                c = mutation(c, bounds, r_mut)
                children.append(c)
        
        for i in range(len(children)):
            # Duplicate 10% of population as previous best candidate
            if i < int(0.1*n_pop):
                children[i] = best

        
        #replace population
        pop = children
    
    return [best, best_eval]





