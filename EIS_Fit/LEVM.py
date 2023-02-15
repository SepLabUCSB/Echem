import sys
import os
import subprocess
import numpy as np
import math
'''

Main function to call for LEVM fitting:

LEVM_fit(freqs, Z, guess, circuit, free_params):
    returns fits: dict of best-fit parameters

'''


def assign_params(circuit, guess, free):
        # Write initial guesses to select parameters
        # Check LEVM manual pp. 114-150 to choose an appropriate
        # function ('A'-'F' most appropriate for solution phase 
        # echem) and map circuit elements to the correct parameter (0-40) 
    
    if circuit == 'Randles_adsorption':
        d = {
             'Rs':   (1, guess['Rs'], free['Rs']),
             'Rct':  (4, guess['Rct'], free['Rct']),
             'Cdl':  (3, guess['Cdl'], free['Cdl']),
             'Cad':  (7, guess['Cad'], free['Cad']),
             'phi':  (9, guess['phi'], free['phi']),
             '_Cad': (10, 2, 0),
             'func': 'C'
             }
        
    elif circuit == 'Randles_uelec':
        d = {
            'R1'  : (1, guess['R1'], free['R1']),
            'R2'  : (4, guess['R2'], free['R2']),
            'R3'  : (21, guess['R3'], free['R3']),
            'Q1'  : (3, guess['Q1'], free['Q1']),
            'Q2'  : (7, guess['Q2'], free['Q2']),
            'n2'  : (9, guess['n2'], free['n2']),
            '_Q2' : (10, 2, 0),
            'func': 'C'            
            }
    
    elif circuit == 'RRC':
        d = {
            'R1'  : (1, guess['R1'], free['R1']),
            'R2'  : (4, guess['R2'], free['R2']),
            'C'   : (3, guess['C'], free['Q1']),
            'func': 'C'
            }
    
    elif circuit == 'RRQ':
        d = {
            'R1'  : (1, guess['R1'], free['R1']),
            'R2'  : (2, guess['R2'], free['R2']),
            'Q1'  : (7, guess['Q1'], free['Q1']),
            'n1'  : (9, guess['n1'], free['n1']),
            '_Q1' : (10, 2, 0),
            'func': 'C'
            }
    
    else:
        print('Circuit not recognized by LEVM.py')
        raise ValueError
        
    return d
        
        


##############################################################################
#####                      LEVM INPUT PARAMS                             #####
##############################################################################

# Input fields defined in LEVM Manual, p. 69
# !! PRESERVE FIELD LENGTHS WHEN CHANGING VALUES !!
inputs = {
    'IOPT':'   0',
    'DINP':'  Z',
    'DFIT':'Z',
    'PINP':'R',
    'PFIT':'R',
    'FREEQ':'F',
    'NEG':' ',
    'FUN':'C',
    'CELCAP':' .0000D+00',
    'DATTYP':'C',
    'IPAR':'         0',
    'ROE':' .0000D+00',
    'IFP':'     7',
    'IRE':'   -11',

    'M':'   27', # M is number of frequencies, automatically determined from input data
    'N':'   40',
    'MAXFEV':'   99',
    'NPRINT':'    0',
    'IRCH':'    3',
    'MODE':'    0',
    'ICP':'    0',
    'IPRINT':'    1',
    'IGACC':'    1',
    'ATEMP':' .0000D+00'
    }




##############################################################################
#####                      GENERAL DATA PROCESSING                       #####
##############################################################################


def float_to_string(n, decimals = 8):
    # Convert scientific notation to Fortran D notation
    if type(n) == str:
        return n
    
    if decimals == 8:
        a = "{0:.7E}".format(n)
    if decimals == 13:
        a = "{0:.12E}".format(n)
    
    digits, exponent = a.split('E')
    
    if digits[0] == '-':
        l = digits[1]
        digits = digits[3:]
        digits = '-.' + l + digits
        
    else:
        l = digits[0]
        digits = digits[2:]
        digits = '.' + l + digits
    
    exponent = int(exponent)
    
    
    if int(exponent) < 0 and int(exponent) > -11:
        s = '%sD-0%s'%(digits, str(abs(exponent)-1))
    
    if int(exponent) < 0 and int(exponent) <= -11:
        s = '%sD-%s'%(digits, str(abs(exponent)-1))
        
    if int(exponent) >= 0 and int(exponent) < 9:
        s = '%sD+0%s'%(digits, str(exponent+1))
    
    if int(exponent) >= 0 and int(exponent) >= 9:
        s = '%sD+%s'%(digits, str(exponent+1))
        
    if digits == '.00000000' and int(exponent) == 0:
        s = '.00000000D+00'
    
    return s



def string_to_float(s):
    # FORTRAN D to scientific notation
    
    digits = float(s[:9])
    sign = s[-3]
    exp = float(s[-2:])
    
    if sign == '+':
        n = digits * 10**(exp)
    
    if sign == '-':
        n = digits * 10**(-1*exp)
    
    return n



def params_to_LEVM_format(params):   
    
    binary_line = [0 for _ in range(40)]
    
    p = {i:0 for i in range(1,41)}
    for key, tup in params.items():
        if key != 'func':
            i, guess, free = tup
            
            p[i] = guess
            binary_line[i-1] = free
            
    binary_line = ''.join(str(j) for j in binary_line)
    function = params['func']
    
    return p, binary_line, function
    



##############################################################################
#####                      INPUT FILE CREATION                           #####
##############################################################################

def write_comment_line(file, comment):
    '''
    Write comment line. 1st line of INFL, 80 char max
    '''
    if len(comment) > 80:
        raise ValueError('Comment is too long!')
        sys.exit()
    
    with open(file, 'w') as f:
        f.write(comment + '\n')
    
    

def write_input_params(file):
    '''
    Write lines 2 and 3 containing fit settings.
    
    Settings defined in global inputs dict
    '''
    
    global inputs
    
    line2 = ''
    line3 = ''
    
    line2_order = ['IOPT', 'DINP', 'DFIT', 'PINP', 'PFIT',
               'FREEQ', 'NEG', 'FUN', 'CELCAP', 'DATTYP',
               'IPAR', 'ROE', 'IFP', 'IRE']


    line3_order = ['M', 'N', 'MAXFEV', 'NPRINT', 'IRCH',
               'MODE', 'ICP', 'IPRINT', 'IGACC', 'ATEMP']

    for key in line2_order:
        line2 = line2 + inputs[key]
    
    for key in line3_order:
        line3 = line3 + inputs[key]
    
    with open(file, 'a') as f:
        f.write(line2 + '\n')
        f.write(line3 + '\n')
     


def write_initial_params(file, p):
    '''
    Write initial guesses to file (lines 4-11)
    '''
    param_lines = {}
    param_lines[1] = '  %s  %s  %s  %s  %s'%(p[1], p[2], p[3], p[4], p[5])
    param_lines[2] = '  %s  %s  %s  %s  %s'%(p[6], p[7], p[8], p[9], p[10])
    param_lines[3] = '  %s  %s  %s  %s  %s'%(p[11], p[12], p[13], p[14], p[15])
    param_lines[4] = '  %s  %s  %s  %s  %s'%(p[16], p[17], p[18], p[19], p[20])
    param_lines[5] = '  %s  %s  %s  %s  %s'%(p[21], p[22], p[23], p[24], p[25])
    param_lines[6] = '  %s  %s  %s  %s  %s'%(p[26], p[27], p[28], p[29], p[30])
    param_lines[7] = '  %s  %s  %s  %s  %s'%(p[31], p[32], p[33], p[34], p[35])
    param_lines[8] = '  %s  %s  %s  %s  %s'%(p[36], p[37], p[38], p[39], p[40])
    
    with open(file, 'a') as f:
        for key, line in param_lines.items():
            f.write(line + '\n')



def write_binary_line(file, binary_line):
    '''
    Line 12 is a 40 character line of binary
    
    If character i == 1, p[i] is free during the fit
    If character i == 0, p[i] is fixed during the fit
    
    '''
    
    with open(file, 'a') as f:
        f.write(binary_line + '\n')



def write_Z_data(file, freqs, Z):
    '''
    Writes lines 13-n, containing all impedance data

    Parameters
    ----------
    file : String.
        Output file. Generally 'INFL'
    freqs : array-like.
        Array of frequencies
    Z : array-like
        Array of impedance data in format (re + 1j*im)

    '''
    freqs = np.asarray(freqs)
    re = np.real(Z)
    im = np.imag(Z)
    
    data_lines = {}

    for i in range(len(freqs)):
        index_val = str(i+1).rjust(5, ' ')
        freq_val = float_to_string(freqs[i], 13).rjust(25, ' ')
        re_val = float_to_string(re[i], 13).rjust(25, ' ')
        im_val = float_to_string(im[i], 13).rjust(25, ' ')
        
        data_lines[i+1] = index_val + freq_val + re_val + im_val
    
    with open(file, 'a') as f:
        for i in data_lines:
            f.write(data_lines[i] + '\n')
    


def write_input_file(file, freqs, Z, params, comment=' '):
    '''
    Parameters
    ----------
    file : String
        Target file to write. Should generally be 'INFL'
    freqs : Array-like
        Array of frequencies.
    Z : Array-like
        Array of (re - 1j*im) impedances.
    params: dict
        dict of {['param'] : (circuit_index, guess, free)}
    comment : String, optional
        Comment line to include on line 1.
        Max 80 characters.

    '''
    
    global inputs
    
    
    p, binary_line, function = params_to_LEVM_format(params)
        
    inputs['FUN'] = function    
    inputs['M'] = str(len(freqs)).rjust(5, ' ')
    
    for i in p:
        p[i] = float_to_string(p[i], 8)
    
    write_comment_line(file, comment)
    write_input_params(file)
    write_initial_params(file, p)
    write_binary_line(file, binary_line)
    write_Z_data(file, freqs, Z)




##############################################################################
#####                      RUN LEVM                                      #####
##############################################################################

def run_LEVM(timeout):
    '''
    Run LEVM using subproccess.run()
    '''
    LEVM_path = os.getcwd() + '\LEVM.EXE'
    try:
        subprocess.run([], executable=LEVM_path, timeout=timeout)
        return 0
    except subprocess.TimeoutExpired:
        print('LEVM.exe timed out')
        return 1


def extract_params(file, params):
    '''
    Function to extract fit parameters from OUTIN
    
    Parameters
    ----------
    file: string
        Output file to read. Should always be 'OUTIN'
    circuit: string
        Name of equivalent circuit used in fit
        (same as in circuits.py)


    Returns
    -------
    d: dict 
        Best-fit parameters, converted back to 
        EIS_fit.py names
    
    ''' 
    with open(file) as f:
        for lineno, line in enumerate(f):
            if lineno == 11:
                b = line
            
    
    p = {}
    for i in range(len(b)):
        if b[i] == '1':
            p[i+1] = 1
            
    with open(file) as f:
        m = 1
        for lineno, line in enumerate(f):
            if lineno > 2 and lineno < 11:
                for element in line.split():
                    p[m] = element
                    m += 1
    
    
    d = {}
    for key, tup in params.items():
        if key == 'func':
            continue
        if key.startswith('_'):
            # Fixed distributed element parameter assignment
            continue
        else:
            i, _, _ = tup
            d[key] = string_to_float(p[i])
          
    
    return d



def LEVM_fit(freqs, Z, guess, circuit, free_params,
             timeout = 2, comment = ' '):
    '''
    Main function to call to perform LEVM fit
    
    Parameters
    ----------
    freqs : Array-like
        List of frequencies.
        
    Z : Array-like
        List of (re - 1j*im) impedances.
        
    guess : Dict
        Initial guesses for fit parameters.


    comment : String, optional
        Comment line to include on line 1.
        Max 80 characters.
        
    Returns
    ---------
    fits : Dict
            Fitted parameters

    '''
    # Determine location of LEVM.py
    path = os.path.realpath(__file__)[:-7]
    os.chdir(path)
    # LEVM should always be in a subfolder called /LEVM/
    # Executable is /LEVM/LEVM.exe
    LEVM_dir = os.getcwd() + '//LEVM/'
    os.chdir(LEVM_dir)
    
    params = assign_params(circuit, guess, free_params)
    
    write_input_file('INFL', freqs, Z, params, comment)
    timedout = run_LEVM(timeout = timeout)
    if timedout == 1:
        return 0
    fits = extract_params('OUTIN', params)
    # Return to LEVM.py directory
    os.chdir(path)
    
    return fits












#%% Testing


## Test fit
# freqs_test = np.array([  100.,   110.,   170.,   250.,   370.,   520.,   750.,  1000.,
#         1300.,  2100.,  3100.,  4400.,  6400.,  9000., 11000., 13000.,
#         17000.])

# Z_test = np.array([18735629.81418155-10792324.49845394j,
#         18106487.57152111-10452511.00087445j,
#         15431234.96575839 -9042780.63944462j,
#         13554600.74026523 -7818661.38105878j,
#         11934174.40651608 -6816407.65898704j,
#         10754284.48652496 -6052968.79942763j,
#         9705710.48484899 -5422119.30431142j,
#         8862943.24254679 -5106026.52578316j,
#         8243324.67947413 -4893743.68317177j,
#         6974518.60516099 -4769951.75226423j,
#         5788202.41140374 -4764130.03690888j,
#         4573546.32794566 -4676062.33443935j,
#         3191922.87380645 -4336168.76449797j,
#         2067134.89945704 -3704618.20436138j,
#         1562232.46402468 -3226321.08913762j,
#         1253126.61199721 -2787895.35068144j,
#         1014951.57014762 -2072276.43707479j])

# d_test = {'R1': 265151.38227128744,
#       'R2': 8488042.354616795,
#       'R3': 41022106.03535404,
#       'Q1': 3.6273344688041505e-12,
#       'n1': 1.0,
#       'Q2': 6.424844265974332e-10,
#       'n2': 0.6790376343555156}


# d_test = {'R1': 2.3911E5,
#       'R2': 9.7956E6,
#       'R3': 3.16791E7,
#       'Q1': 4.07581E-12,
#       'n1': 1.0,
#       'Q2': 4.93661E-10,
#       'n2': 0.704902}



if __name__ == '__main__':
    
    # file = r'C:/Users/BRoehrich/Desktop/EIS-EAB data/2022-09-14/vanco -340mV/0006s.txt'

    # f, re, im = np.loadtxt(file, skiprows=1, unpack=True)
    # Z = re + 1j*im
    
    # d_test = {'Rs': 8.7437e+002,
    #           'Rct': 1.94050e+004,
    #           'Cad': 9.99434e-007,
    #           'phi':8.40000e-001,
    #           'Cdl':3.69744e-007}
    
    # free_params = {'Rs': 0,
    #           'Rct': 1,
    #           'Cad': 1,
    #           'phi':0,
    #           'Cdl':1}
    
    # fits = LEVM_fit(f, Z, d_test, 'Randles_adsorption', free_params)
    # print(fits)
    # a = os.path.realpath(__file__)
    
    # folder = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-09-14\vanco -340mV'
    folder = r'C:\Users\BRoehrich\Desktop\EIS-EAB data\2022-11-16\first 30 min -340mV'
    params = []
    kets = []
    # d = {'Rs': 8.7437e+002,
    #         'Rct': 1.94050e+004,
    #         'Cad': 9.99434e-007,
    #         'phi':8.40000e-001,
    #         'Cdl':3.69744e-007}
    
    d = {'Rs': 543,
            'Rct': 3.44050e+004,
            'Cad': 4.79434e-007,
            'phi':8.40000e-001,
            'Cdl':2.25e-007}   
    
    free_params = {'Rs': 1,
            'Rct': 1,
            'Cad': 1,
            'phi':0,
            'Cdl':1}
    
    
    i = 0
    for file in os.listdir(folder):
        if file.endswith('s.txt'):
            while i < 450:
                if file.endswith('currents.txt'):
                    continue
                if file.endswith('s.txt'):
                    f, re, im = np.loadtxt(os.path.join(folder,file), skiprows=1, unpack=True)
                    Z = re + 1j*im
                    
                    if i == 0:
                        d = d
                    else:
                        d = params[i-1]
                    # d['Rs'] = re[-1]
                    
                    fits = LEVM_fit(f, Z, d, 'Randles_adsorption', free_params)
                    print(fits)
                    
                    params.append(fits)
                    ket = 1/(2*fits['Rct']*fits['Cad'])
                    kets.append(ket)
                    
                    i += 1
                                    
                            
    

