import numpy as np
import os
import sys
import subprocess

# Set path to dStar directory to import reader for reading history.data files
home_directory = os.getenv('HOME')
dStar_directory = home_directory + '/dStar/'
sys.path.append(dStar_directory + 'tools/')

import reader

def log_probability(parameter_values, Q_heating_inner,
                    parameter_minimum_values=np.array([8.35e7, 0.0, 0.0]), 
                    parameter_maximum_values=np.array([9.35e8, 20.0, 20.0])):
    '''
    Calculates the log of the probability of a model with dStar model 
    parameters for use in an MCMC.
    
    The log of the probability is based on the chi^2 error assuming a Gaussian
    probability distribution: probability = exp(-chi2/2). I ignore leading
    coefficients since the MCMC only uses ratios of probabilities so they would
    cancel anyway.
    
    Uses uniform priors for all model parameters with given minimum values and 
    maximum values.
    
    Parameters
    ----------
    
    parameter_values: array-like'
        Input array of parameters for dStar model.
        Must modify the src/run.f file to accomodate wrappers for the chosen 
        parameters.
        Parameters here are: [core_temperature, Qimp, Q_heating_shallow]
        
    parameter_minimum_values: array-like
        Minimum values for each dStar model parameter.
        If any parameter is less than its corresponding minimum value, this 
        function will return -np.inf (0 probability).
        
    parameter_maximum_values: array-like
        Maximum values for each dStar model parameter.
        If any parameter is greater than its corresponding minimum value, this 
        function will return -np.inf (0 probability).
        
    Returns
    -------
    
    log_probability: float
        The log of the probability for use in the MCMC sample selection 
        process.
    
    '''
    # Check for any parameters less than the minimum values or greater than the 
    # maximum values.
    # If any parameters are outside the allowed range, then return -np.inf.
    if np.any(parameter_values < parameter_minimum_values) or \
        np.any(parameter_values > parameter_maximum_values):
        return -np.inf
    # If paramaters are within the allowed range, run the models and calculate 
    # the log of the probability.
    else:
        # Run the model and save the output chi2
        # Convert all parameter values to strings so they can be used as input 
        # in the subprocess command.
        core_temperature = str(parameter_values[0])
        Qimp = str(parameter_values[1])
        Q_heating_shallow = str(parameter_values[2])
        # Try to run the dStar model which calculates and outputs the chi2.
        # If the model can't run or if the chi2 is too high, it won't produce
        # output that can be converted to a float and will return an error.
        # In this case just return -np.inf (probability = 0).
        try:
            chi2 = float(subprocess.run(['./run_dStar', '-T', core_temperature, 
                                         '-Q', Qimp, '-H', Q_heating_shallow], 
                                        capture_output=True, text=True).stdout)
            return -chi2/2
        except:
            return -np.inf
        
def load_LOGS_data(LOGS_directory):
    
    f = reader.dStarReader('{}/history.data'.format(LOGS_directory))

    # time after outburts in days
    t = f.data.time
    # log(L) [ergs]
    lg_L = f.data.lg_Lsurf
    # log(T_eff) [K]
    lg_Teff = f.data.lg_Teff
    
    # constants in cgs
    c = 3e10
    G = 6.67e-8
    M_sun = 1.989e33

    M = M_sun*f.header.total_mass
    R = f.header.total_radius*1e5
    redshift_factor = (1 - 2*G*M/(R*c**2))**(-.5)
    
    Teff_inf = 10**(lg_Teff)/redshift_factor
    
    return t, lg_L, lg_Teff, Teff_inf

def write_inlist(parameter_values, parameter_names, base_inlist_name='inlist_base', new_inlist_name='inlist', LOGS_directory='LOGS'):

    base_inlist = open(base_inlist_name, 'r')
    new_inlist = open(new_inlist_name, 'w')

    base_inlist.seek(0)

    for line in base_inlist.readlines():

        parameter_line = False
        for i, parameter_name in enumerate(parameter_names):
            if line[:7+len(parameter_name)] == '    {} = '.format(parameter_name):
                if parameter_name == 'core_temperature':
                    new_inlist.write('    {} = {:.4e}\n'.format(parameter_name, parameter_values[i]))
                    parameter_line = True
                    break
                else:
                    new_inlist.write('    {} = {:.4f}\n'.format(parameter_name, parameter_values[i]))
                    parameter_line = True
                    break

        # If we just modified a parameter line, move to the next line.
        if parameter_line == True:
            continue
        
        # When we reach the output_directory assignment line, reassign to the new output_directory.
        elif line[:23] == '    output_directory = ':
            new_inlist.write("    output_directory = '{}'\n".format(LOGS_directory))

        # If it is not one of the lines to change, then copy it to each new inlist.
        else:
            new_inlist.write(line)
            
    base_inlist.close()
    new_inlist.close()