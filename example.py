import pycali
import matplotlib.pyplot as plt 
import numpy as np

#######################################################
# setup configurations, there are two ways:
# 1) load from a param file
#    cfg = pycali.Config("param.txt")
# 2) direct call setup()
# 
cfg = pycali.Config()

# except for the argument "fcont", the rest arguments are optional.
# e.g.,  cfg.setup(fcont="data/ngc5548_cont.txt")
#
cfg.setup(
          fcont="data/ngc5548_cont.txt",     # fcont is a string 
          fline=["data/ngc5548_line.txt"],   # fline is a list, include multiple lines
          nmcmc=10000, ptol=0.1,
          scale_range_low=0.5, scale_range_up=2.0,
          shift_range_low=-1.0, shift_range_up=1.0,
          syserr_range_low=0.0, syserr_range_up=0.2,
          errscale_range_low=0.5, errscale_range_up=2.0,
          sigma_range_low=1.0e-4, sigma_range_up=1.0,
          tau_range_low=1.0, tau_range_up=1.0e4,
          fixed_scale=False, fixed_shift=False,
          fixed_syserr=True, fixed_error_scale=True,
          fixed_codes=[], # fixed_codes is a list to specify the codes that need not to intercalibrate
                          # e.g., [1, 3], will fix 1st and 3rd codes
          fixed_scalecodes=[], # fixed_scalecodes is a list to specify the codes that need to fix scale (to 1)
                          # e.g., [1, 3], will fix scale of 1st and 3rd codes
          flag_norm=True, # whether do normalization before intercalibrating
          )
cfg.print_cfg()

######################################################
# do intercalibration
#
cali = pycali.Cali(cfg)  # create an instance
cali.mcmc()              # do mcmc
cali.get_best_params()   # calculate the best parameters
cali.output()            # print output
cali.recon()             # do reconstruction

# plot results to PyCALI_results.pdf
# set smooth=True for smooth histogram plots
pycali.plot_results(cfg, smooth=False)

# a simple plot 
pycali.simple_plot(cfg)
