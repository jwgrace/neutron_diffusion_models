import numpy as np
import os
import sys
import emcee
import corner
import matplotlib.pyplot as plt

# Set path to dStar directory to import reader for reading history.data files
home_directory = os.getenv('HOME')
dStar_directory = home_directory + '/dStar/'
sys.path.append(dStar_directory + 'tools/')

import reader as dStar_reader

h5_filename = 'test.h5'
parameter_names = ['Core Temperature (K)', 'Impurity Parameter', 'Shallow Heating (MeV)']
corner_plot_title = 'KS 1701-260 (Traditional)'
corner_plot_name = 'corner_plot-traditional.png'

ndim = len(parameter_names)

reader = emcee.backends.HDFBackend(h5_filename)

tau = reader.get_autocorr_time()
burnin = int(2 * np.max(tau))
#thin = int(0.5 * np.min(tau))
samples = reader.get_chain(discard=burnin, flat=True)
log_prob_samples = reader.get_log_prob(discard=burnin, flat=True)
log_prior_samples = reader.get_blobs(discard=burnin, flat=True)

fig = corner.corner(samples, labels=parameter_names)
fig.patch.set_facecolor('white')
fig.set_dpi(200)

fig.suptitle(corner_plot_title, fontsize=18)
fig.savefig(corner_plot_name)

# Arrays for the best-fit parameters from the MCMC as well as the upper and lower errors.
parameter_values_mcmc = np.zeros_like(parameter_values_initial)
err_low = np.zeros_like(parameter_values_mcmc)
err_high = np.zeros_like(parameter_values_mcmc)

for i in range(ndim):
    mcmc = np.percentile(samples[:, i], [16, 50, 84])
    parameter_values_mcmc[i] = mcmc[1]
    err_low[i] = mcmc[1] - mcmc[0]
    err_high[i] = mcmc[2] - mcmc[1]
    
    if parameter_names[i] == 'Core Temperature (K)':
        print('{}\t {0:.3e} -{:.3e} + {:.3e}'.format(parameter_names[i], mcmc[i], err_low[i], err_high[i]))
    else:
        print('{}\t {0:.3f} -{:.3f} + {:.3f}'.format(parameter_names[i], mcmc[i], err_low[i], err_high[i]))
    
chi2 = -2*log_probability(parameter_values_mcmc)
print('chi2 =', chi2)

degrees_of_freedom = len(t_obs) - ndim
reduced_chi2 = chi2/degrees_of_freedom
print(reduced_chi2)