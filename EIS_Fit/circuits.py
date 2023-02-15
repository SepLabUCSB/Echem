import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('Z:/Projects/Brian/scientific.mplstyle')


##############################################################################
#####                        DEFINE CIRCUITS                             #####
##############################################################################

def list_circuits():
    #
    # !! Add circuits here when they are defined below !!
    #
    
    l = ['CPE', 'RC',  'RRC', 'RRQ', 'Randles', 'Randles_CPE', 
                    'Randles_uelec', 'Randles_adsorption']
    funcs = [CPE, RC, RRC, RRQ, Randles, Randles_CPE, Randles_uelec,
             Randles_adsorption]
    return l, funcs


def CPE(f, params):
    '''
    Params:
        Q: CPE value
        n: CPE exponent (0 < n < 1)
    '''
    w = 2*np.pi*f
    Q = params['Q']
    n = params['n']
    return 1/(Q*(w*1j)**n)


def RC(f, params):
    '''
    Params:
        R: Resistance
        C: Capacitance
    '''
    w = 2*np.pi*f
    R = params['R']
    C = params['C']
    
    Z_C = 1/(1j*w*C)
    return (R*Z_C)/(Z_C + R)

def RQ(f, params):
    '''
    Params:
        R: Resistance
        Q: CPE value
        n: CPE exponent (0 < n < 1)
    '''
    w = 2*np.pi*f
    R = params['R']
    Q = params['Q']
    n = params['n']
    
    Z_C = 1/(Q*(w*1j)**n)
    return (R*Z_C)/(Z_C + R)


def RRC(f, params):
    '''
    Params:
        R1: Series resistance
        R2: Parallel resistance
        C: Capacitance
    '''
    w = 2*np.pi*f
    R1 = params['R1']
    R2 = params['R2']
    C = params['C']
    
    Z_C = 1/(1j*w*C)
    
    return R1 + (R2*Z_C)/(Z_C + R2)


def RRQ(f, params):
    '''
    Params:
        R1: Series resistance
        R2: Parallel resistance
        Q1: Capacitance CPE
        n1: CPE phase
    '''
    w = 2*np.pi*f
    R1 = params['R1']
    R2 = params['R2']
    Q1 = params['Q1']
    n1 = params['n1']
    
    Z_C = 1/(Q1*(w*1j)**n1)
    
    return R1 + (R2*Z_C)/(Z_C + R2)


def Randles(f, params):
    '''
    Params:
        R1: Series resistance
        R2: CT resistance
        C: dl capacitance
        sigma: Warburg sigma

    '''
    # R-(C, RW)
    w = 2*np.pi*f
    R1 = params['R1']
    R2 = params['R2']
    C = params['C']
    sigma = params['sigma']
    
    Z_C = 1/(1j*w*C)
    Z_W = sigma*(w**(-0.5))-1j*sigma*(w**(-0.5))
    
    return R1 + 1/(1/Z_C + 1/(R2+Z_W)) 


def Randles_CPE(f, params):
    '''
    Params:
        R1: Series resistance
        R2: CT resistance
        Q: dl capacitance
        n: dl exponent (0 < n < 1)
        sigma: Warburg sigma

    '''
    w = 2*np.pi*f
    R1 = params['R1']
    R2 = params['R2']
    Q = params['Q']
    n = params['n']
    sigma = params['sigma']
    
    Z_Q = 1/(Q*(w*1j)**n)
    Z_W = sigma*(w**(-0.5))-1j*sigma*(w**(-0.5))
    
    return R1 + 1/(1/Z_Q + 1/(R2+Z_W))     


def Randles_uelec(f, params):
    '''
    Params:
        R1: Series resistance
        R2: CT resistance
        R3: Diff. resistance
        Q1: dl capacitance
        n1: dl exponent (0 < n < 1)
        Q2: Warburg CPE
        n2: Warburg exponent (0 < n < 1)
    '''

    R1 = params['R1']
    R2 = params['R2']
    R3 = params['R3']
    Q1 = params['Q1']
    n1 = params['n1']
    Q2 = params['Q2']
    n2 = params['n2']
        
    
    Z_d = RQ(f,{'R':R3, 'Q':Q2, 'n':n2})
    C_dl = CPE(f, {'Q':Q1, 'n':n1})
    
    Z_top = C_dl*(R2 + Z_d)
    Z_bottom = C_dl + R2 + Z_d
    
    Z_tot = R1 + Z_top/Z_bottom

    return Z_tot



def Randles_adsorption(f, params):
    '''
    Params:
        Rs: Series resistance
        Rct: CT resistance
        Cdl: dl capacitance
        ndl: dl exponent (0 < n < 1)
        Cad: Adsorption CPE
        nad: Adsorption exponent (0 < n < 1)
    '''
    
    R1 = params['Rs']
    R2 = params['Rct']
    Q1 = params['Cdl']
    n1 = params['ndl']
    Q2 = params['Cad']
    n2 = params['nad']
    
    Ca = CPE(f, {'Q':Q2, 'n':n2})
    Cdl = CPE(f, {'Q':Q1, 'n':n1})
    
    Z = R1 + 1/(1/Cdl + 1/(R2+Ca))
    
    return Z
    


############################################
##      Add new circuits above here!      ##
############################################

def chosen_circuit(circuit, w, params):
    # Helper function for iteratively plotting fits
    circuit_list, circuit_funcs = list_circuits()
    
    circuitfunc = circuit_funcs[circuit_list.index(circuit)]
    
    return circuitfunc(w, params)







##############################################################################
#####              DEFINE ERROR FUNCTIONS TO MINIMIZE                    #####
############################################################################## 

def leastsq_errorfunc(w, Z, params, circuit):
    '''

    Parameters
    ----------
    w : array of frequencies
    Z : array of (re + im)
    params : dict of fit parameters.
    circuit : string representing circuit

    Returns
    -------
    S : Weighted sum of squares.

    '''
    
    circuit_list, circuit_funcs = list_circuits()
    
    
    if not circuit in circuit_list:
        print('Circuit not recognized. Allowed circuits: ')
        for c in circuit_list:
            print(c)
        raise ValueError('Circuit not listed in circuits.py')
    
    circuitfunc = circuit_funcs[circuit_list.index(circuit)]
    
    
    
    Z = np.asarray(Z)

    re = np.real(Z)
    im = np.imag(Z)
    
    
    Z_fit = circuitfunc(w, params)
    
    re_fit = np.real(Z_fit)
    im_fit = np.imag(Z_fit) 
    
    error = np.array([(re-re_fit)**2, (im-im_fit)**2])
    
    
    # Modulus weighting (proportional to 1/Z)
    weight = np.array([1/((re_fit**2 + im_fit**2)**(1/2)), 1/((re_fit**2 + im_fit**2)**(1/2))])
    
    # weight = np.ones(len(re_fit))
    
    S = np.sum(weight*error)
            
    return S



def calc_chi_squared(w, Z, params, circuit):
    '''
    Copy of leastsq_errorfunc that returns chi squared instead
    
    Unweighted so not preferred for fitting
    
    Parameters
    ----------
    w : array of frequencies
    Z : array of (re + im)
    params : dict of fit parameters.
    circuit : string representing circuit

    Returns
    -------
    S : Weighted sum of squares.

    '''
    
    circuit_list, circuit_funcs = list_circuits()
    
    
    if not circuit in circuit_list:
        print('Circuit not recognized. Allowed circuits: ')
        for c in circuit_list:
            print(c)
        raise ValueError('Circuit not listed in circuits.py')
    
    circuitfunc = circuit_funcs[circuit_list.index(circuit)]
    
    
    
    Z = np.asarray(Z)

    re = np.real(Z)
    im = np.imag(Z)
    
    
    Z_fit = circuitfunc(w, params)
    
    re_fit = np.real(Z_fit)
    im_fit = np.imag(Z_fit) 
    
    error = np.array([(re-re_fit)**2, np.abs((im-im_fit)**2)])
    
    chi_squared = np.sum(error/np.array([np.abs(Z)]))
    
    chi_squared = chi_squared/len(w)
        
    return chi_squared
        


def leastsq_errorfunc_array(param_array, w, Z, circuit, param_names):
    '''
    Array-based version of leastsq_errorfunc for
    compatibility with scipy.minimize

    Parameters
    ----------
    param_array : array
        Array of equivalent circuit parameter values.
    w : array
        Array of frequencies.
    Z : array
        Array of impedance data, re + 1j*im.
    circuit : string
        String description of circuit from circuits.py
    param_names : list
        list of parameter names to convert back to params dict


    Returns
    -------
    S : float
        Weighted sum of squares.

    '''
    
    circuit_list, circuit_funcs = list_circuits()
    
    
    if not circuit in circuit_list:
        print('Circuit not recognized. Allowed circuits: ')
        for c in circuit_list:
            print(c)
        raise ValueError('Circuit not listed in circuits.py')
    
    
    circuitfunc = circuit_funcs[circuit_list.index(circuit)]
    
    params = dict()
    for i in range(len(param_array)):
        name = param_names[i]
        val = param_array[i]
        params[name] = val
    
    Z = np.asarray(Z)

    re = np.real(Z)
    im = np.imag(Z)
    
    
    Z_fit = circuitfunc(w, params)
    
    re_fit = np.real(Z_fit)
    im_fit = np.imag(Z_fit) 
    
    error = np.array([(re-re_fit)**2, (im-im_fit)**2])
    
    
    # Modulus weighting (proportional to 1/Z)
    weight = np.array([1/((re_fit**2 + im_fit**2)**(1/2)), 1/((re_fit**2 + im_fit**2)**(1/2))])
    
    # weight = np.ones(len(re_fit))
    
    S = np.sum(weight*error)
        
    return S





