import sys
import os
import subprocess
import numpy as np
import math
'''

Main function to call for LEVM fitting:

LEVM_fit(freqs, Z, d, circuit):
    returns fits: dict of best-fit parameters

'''


##############################################################################
#####                      LEVM INPUT PARAMS                             #####
##############################################################################

# Input fields defined in LEVM Manual, p. 69
# !! PRESERVE FIELD LENGTHS WHEN CHANGING PARAMETERS !!
inputs = {
    'IOPT':'   0',
    'DINP':'  Z',
    'DFIT':'Z',
    'PINP':'Z',
    'PFIT':'R',
    'FREEQ':'F',
    'NEG':' ',
    'FUN':'C',
    'CELCAP':' .0000D+00',
    'DATTYP':'C',
    'IPAR':'         0',
    'ROE':' .0000D+00',
    'IFP':'     1',
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



def params_to_LEVM_format(d, circuit):
    '''
    Map initial guesses to correct parameters[0-40]

    Parameters
    ----------
    d : Dict
        Initial parameter guesses (same as for circuits.py)
    circuit : String
        Name of equivalent circuit (same as for circuits.py)


    Returns
    -------
    p: dict, 
        p[0] to p[40]
    function: string,
        Built-in LEVM equivalent circuit to use

    '''
    
    def _round(x):
        OOM = math.floor(math.log10(x))
        # return round(x/(10**OOM))*10**OOM
        # return 10**OOM
        return x
    
    
    if circuit == 'Randles_uelec':
        # Write initial guesses to select parameters
        # Check LEVM manual pp. 114-150 to choose an appropriate
        # function ('A'-'F' most appropriate for solution phase 
        # echem) and map circuit elements to the correct parameter (0-40)
        #
        # Circuits must also be added to extract_params() below!!
        p = {}
        p[1] = d['R1']
        p[4] = d['R2']
        p[21] = d['R3']
        p[9] = d['n2']
        p[3] = d['Q1']
        p[7] = d['Q2']
        p[10] = 2
        function = 'C'
    
    
    elif circuit == 'Randles_adsorption':
        p = {}
        p[1] = d['R1']  # Rs
        p[24] = d['R2']  # Rct
        p[12] = d['Q1']  # Cdl
        p[14] = d['n1']  # Cdl phase
        p[15] = 2        # Assign CPE for Cdl
        p[17] = d['Q2']  # Cad
        p[19] = d['n2']  # Cad phase
        p[20] = 2       # Assign CPE for Cad
        function = 'C'
    
    
    elif circuit == 'RRC':
        p = {}
        p[1] = d['R1']
        p[4] = d['R2']
        p[3] = d['C']
        function = 'C'
    
        
    elif circuit == 'RRQ':
        p = {}
        p[1] = d['R1']
        p[2] = d['R2']
        p[7] = d['Q1']
        p[9] = d['n1']
        p[10] = 2
        function = 'C'
        
    else:
        print('Circuit not recognized LEVM line 174')
        raise ValueError
    

    
    for i in range(1,41):
        if i not in p:
            p[i] = 0
    
    
    return p, function
    



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



def write_binary_line(file, p):
    '''
    Line 12 is a 40 character line of binary
    
    If character i == 1, p[i] is free during the fit
    If character i == 0, p[i] is fixed during the fit
    
    '''
    
    line = ''
    
    keys = list(p.keys())
    keys.sort()
    for key in keys:
        val = string_to_float(p[key])
        
        if val == 2:
            # CPE assigned using NLEM == 2, fixed value
            line = line + '0'
            
        if val == 1:
            # Used to force CPE to be a capacitor
            # Not a free parameter
            line = line + '0'
        
        elif val != 0:
            line = line + '1'
            
        elif val == 0:
            line = line + '0'
    
    with open(file, 'a') as f:
        f.write(line + '\n')



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
    


def write_input_file(file, d, freqs, Z, circuit, comment=' '):
    '''
    Parameters
    ----------
    file : String
        Target file to write. Should generally be 'INFL'
    d : Dict
        Initial guesses for fit parameters.
    freqs : Array-like
        Array of frequencies.
    Z : Array-like
        Array of (re - 1j*im) impedances.
    comment : String, optional
        Comment line to include on line 1.
        Max 80 characters.

    '''
    
    global inputs
    
    
    p, circuit = params_to_LEVM_format(d, circuit)
        
    inputs['FUN'] = circuit    
    inputs['M'] = str(len(freqs)).rjust(5, ' ')
    
    for i in p:
        p[i] = float_to_string(p[i], 8)
    
    write_comment_line(file, comment)
    write_input_params(file)
    write_initial_params(file, p)
    write_binary_line(file, p)
    write_Z_data(file, freqs, Z)




##############################################################################
#####                      RUN LEVM                                      #####
##############################################################################

def run_LEVM(timeout):
    '''
    Run LEVM using subproccess.run()
    '''
    LEVM_path = os.getcwd() + '\LEVM.EXE'
    subprocess.run([], executable=LEVM_path, timeout=timeout)
           


def extract_params(file, circuit):
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
                    if m in p.keys():
                        p[m] = element
                    m += 1
                    
    if circuit == 'Randles_uelec':
        d = {}
        d['R1'] = string_to_float(p[1])
        d['Q1'] = string_to_float(p[3])
        d['n1'] = 1.0
        d['R2'] = string_to_float(p[4])
        d['Q2'] = string_to_float(p[7])
        d['n2'] = string_to_float(p[9])
        d['R3'] = string_to_float(p[21])
    
    
    if circuit == 'Randles_adsorption':
        d = {}
        d['R1'] = string_to_float(p[1])
        d['R2'] = string_to_float(p[4])
        d['n1'] = 1.0    
        d['Q1'] = string_to_float(p[3])
        d['Q2'] = string_to_float(p[7])
        d['n2'] = string_to_float(p[9])
    
    
    if circuit == 'RRC':
        d = {}
        d['R1'] = string_to_float(p[1])
        d['R2'] = string_to_float(p[4])
        d['C'] = string_to_float(p[3])
        
        
    if circuit == 'RRQ':
        d = {}
        d['R1'] = string_to_float(p[1])
        d['R2'] = string_to_float(p[2])   
        d['Q1'] = string_to_float(p[7])
        d['n1'] = string_to_float(p[9])
        
        
    else:
        print('Circuit not recognized')
        raise ValueError('Circuit not recognized by extract_params()')
    
    return d



def LEVM_fit(freqs, Z, d, circuit, timeout = 2, comment = ' '):
    '''
    Main function to call to perform LEVM fit
    
    Parameters
    ----------
    freqs : Array-like
        List of frequencies.
        
    Z : Array-like
        List of (re - 1j*im) impedances.
        
    d : Dict
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
    
    write_input_file('INFL', d, freqs, Z, circuit, comment)
    run_LEVM(timeout = timeout)
    fits = extract_params('OUTIN', circuit)
    
    # Return to LEVM.py directory
    os.chdir(path)
    
    return fits












#%% Testing


## Test fit
freqs_test = np.array([  100.,   110.,   170.,   250.,   370.,   520.,   750.,  1000.,
        1300.,  2100.,  3100.,  4400.,  6400.,  9000., 11000., 13000.,
        17000.])

Z_test = np.array([18735629.81418155-10792324.49845394j,
        18106487.57152111-10452511.00087445j,
        15431234.96575839 -9042780.63944462j,
        13554600.74026523 -7818661.38105878j,
        11934174.40651608 -6816407.65898704j,
        10754284.48652496 -6052968.79942763j,
        9705710.48484899 -5422119.30431142j,
        8862943.24254679 -5106026.52578316j,
        8243324.67947413 -4893743.68317177j,
        6974518.60516099 -4769951.75226423j,
        5788202.41140374 -4764130.03690888j,
        4573546.32794566 -4676062.33443935j,
        3191922.87380645 -4336168.76449797j,
        2067134.89945704 -3704618.20436138j,
        1562232.46402468 -3226321.08913762j,
        1253126.61199721 -2787895.35068144j,
        1014951.57014762 -2072276.43707479j])

d_test = {'R1': 265151.38227128744,
      'R2': 8488042.354616795,
      'R3': 41022106.03535404,
      'Q1': 3.6273344688041505e-12,
      'n1': 1.0,
      'Q2': 6.424844265974332e-10,
      'n2': 0.6790376343555156}


# d_test = {'R1': 2.3911E5,
#       'R2': 9.7956E6,
#       'R3': 3.16791E7,
#       'Q1': 4.07581E-12,
#       'n1': 1.0,
#       'Q2': 4.93661E-10,
#       'n2': 0.704902}

if __name__ == '__main__':
    fits = LEVM_fit(freqs_test, Z_test, d_test, 'Randles_uelec')
    print(fits)
    # a = os.path.realpath(__file__)
    

