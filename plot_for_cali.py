#
# plots for cali 
# 
# run "python plot_for_cali.py" after "./cali param.txt"
#

import pycali
import matplotlib.pyplot as plt 
import numpy as np

# load configuration form param.txt
cfg = pycali.Config("param.txt")

# print cfg
cfg.print_cfg()

# plot results to PyCALI_results.pdf
pycali.plot_results(cfg)

# a simple plot 
pycali.simple_plot(cfg)