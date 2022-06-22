import numpy as np
import sys
import getopt
import subprocess
import emcee
from multiprocessing import Pool
import ns_modelling as ns

argument_list = sys.argv[1:]
options = 'i:n:w:o:'
long_options = ['Q_heating_inner=', 'max_steps=', 'nwalkers=', 'output_file=']

# Default values for inner heating and MCMC parameters.
# These will be overwritten if values are supplied by command-line options.
Q_heating_inner = 1.5
max_n = 10000
nwalkers = 32
h5_filename = 'test.h5'

try:
    # Parsing argument
    arguments, values = getopt.getopt(argument_list, options, long_options)

    for argument, value in arguments:
        if argument in ('-i', '--Q_heating_inner'):
            try:
                Q_heating_inner = float(value)
            except ValueError:
                print('Q_heating_inner must be a float')
        if argument in ('-n', '--max_steps'):
            try:
                max_n = int(value)
            except ValueError:
                print('max_n must be an integer')
        if argument in ('-w', '--nwalkers'):
            try:
                nwalkers = int(value)
            except ValueError:
                print('nwalkers must be an integer')
        if argument in ('-o', '--output_file'):
            try:
                h5_filename = str(value)
            except ValueError:
                print('h5_filename must be a string')

except getopt.error as err:
    # output error, and return with an error code
    print (str(err))

def log_probability(parameter_values, Q_heating_inner=1.5,
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
                                         '-Q', Qimp, '-s', Q_heating_shallow,
                                         '-i', Q_heating_inner],
                                        capture_output=True, text=True).stdout)
            return -chi2/2
        except:
            return -np.inf

# Initial guesses for parameters: Core Temperature, Qimp, Q_heating_shallow
# Start sampling from fit_lightcurve example in dStar directory
parameter_values_initial = np.array([9.35e7, 4.2, 1.37])

ndim = len(parameter_values_initial)

# 2D array of initial positions of the walkers.
# Center the walkers on the initial guess and randomly vary each parameter by up to 10%.
pos = parameter_values_initial + .1*parameter_values_initial*np.random.randn(nwalkers, ndim)

backend = emcee.backends.HDFBackend(h5_filename)
backend.reset(nwalkers, ndim)

with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, kwargs=(Q_heating_inner=Q_heating_inner), pool=pool, backend=backend)
    #sampler.run_mcmc(pos, max_steps, progress=True)

    # This will be useful to testing convergence
    old_tau = np.inf

    # Now we'll sample for up to max_n steps
    for sample in sampler.sample(pos, iterations=max_n, progress=True):
        # Only check convergence every 100 steps
        if sampler.iteration % 100:
            continue

        # Compute the autocorrelation time so far
        # Using tol=0 means that we'll always get an estimate even
        # if it isn't trustworthy
        tau = sampler.get_autocorr_time(tol=0)

        # Check convergence
        converged = np.all(tau * 100 < sampler.iteration)
        converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
        if converged:
            break
        old_tau = tau
