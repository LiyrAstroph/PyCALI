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
# 
cfg.setup(fcont="data/ngc5548_cont.txt", fline="data/ngc5548_line.txt", 
          nmcmc=10000, pdiff=0.1,
          scale_range_low=0.5, scale_range_up=1.5,
          shift_range_low=-1.0, shift_range_up=1.0,
          sigma_range_low=1.0e-4, sigma_range_up=1.0,
          tau_range_low=1.0, tau_range_up=1.0e4)
cfg.print_cfg()

######################################################
# do intercalibration
#
cali = pycali.Cali(cfg)  # create an instance
cali.mcmc()              # do mcmc
cali.get_best_params()   # calculate the best parameters
cali.align_with_error()  # align the light curves
cali.output()            # print output
cali.recon()             # do reconstruction

######################################################
# now plot
# 
data={}
nax = 1
cont = np.loadtxt(cfg.fcont)
cont_cali = np.loadtxt(cfg.fcont+"_cali")
data["cont"]=[cont, cont_cali]

if cfg.fline:
  nax+=1
  line = np.loadtxt(cfg.fline)
  line_cali = np.loadtxt(cfg.fline+"_cali")
  data["line"] = [line, line_cali]

fig = plt.figure()
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
for i, key in enumerate(data.keys()):
  ax = fig.add_subplot(nax, 1, i+1)
  d = data[key][0]
  dc = data[key][1]
  ax.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], ls='none', marker='o', markersize=4, color=cycle[0], 
              ecolor='darkgrey', markeredgecolor=None, elinewidth=1, label=key)
  ax.errorbar(dc[:, 0], dc[:, 1], yerr=dc[:, 2], ls='none', marker='o', markersize=4, color=cycle[1],
              ecolor='darkgrey', markeredgecolor=None,  elinewidth=1, label=key+"_cali")
  ax.legend()
  ax.set_xlabel("Time")
  ax.set_ylabel("Flux")

plt.show()


