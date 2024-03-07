#
# plots for cali 
# 
# run "python plot_for_cali.py param.txt" after "./cali param.txt"
#
# this file should be placed in the same directory where "data/" exists
#

import sys
import corner
import configparser as cfgpars 
import matplotlib.pyplot as plt 
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from os.path import basename
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

EPS = np.finfo(np.float64).tiny

class Config:
  def __init__(self, fname="param.txt"):
    config = cfgpars.RawConfigParser(delimiters=' ', comment_prefixes='#', inline_comment_prefixes='#', \
      default_section=cfgpars.DEFAULTSECT, empty_lines_in_values=False)
    
    with open(fname) as f:
      file_content = '[dump]\n' + f.read()

    config.read_string(file_content)

    param = config["dump"]
    self.fcont = param["FileCont"]

    if "FileLine" in param:
      self.fline = param["FileLine"].split(",")
      # strip the begining and end space
      for i in range(len(self.fline)):
        self.fline[i] = self.fline[i].strip()
    else:
      self.fline = []
    
    if "NMCMC" in param:
      self.nmcmc = int(param["NMCMC"])
    else:
      self.nmcmc = 10000
    
    if "PTol" in param:
      self.ptol = float(param["PTol"])
    else:
      self.ptol = 0.1
    
    if "FixedScale" in param:
      self.fixed_scale = int(param["FixedScale"])
    else:
      self.fixed_scale = 0
    
    if "FixedShift" in param:
      self.fixed_scale = int(param["FixedShift"])
    else:
      self.fixed_scale = 0

    if "FixedSyserr" in param:
      self.fixed_syserr = int(param["FixedSyserr"])
    else:
      self.fixed_syserr = 1

    if "FixedErrorScale" in param:
      self.fixed_error_scale = int(param["FixedErrorScale"])
    else:
      self.fixed_error_scale = 1

    # parameter range 
    # ===================Simga=============================
    if "SigmaRangeLow" in param:
      self.sigma_range_low = float(param["SigmaRangeLow"])
    else:
      self.sigma_range_low = 1.0e-4
    
    if "SigmaRangeUp" in param:
      self.sigma_range_up = float(param["SigmaRangeUp"])
    else:
      self.sigma_range_up = 1.0
    
    # ===================Tau=============================
    if "TauRangeLow" in param:
      self.tau_range_low = float(param["TauRangeLow"])
    else:
      self.tau_range_low = 1.0
    
    if "TauRangeUp" in param:
      self.tau_range_up = float(param["TauRangeUp"])
    else:
      self.tau_range_up = 1.0e4
    
    # ===================Scale=============================
    if "ScaleRangeLow" in param:
      self.scale_range_low = float(param["ScaleRangeLow"])
    else:
      self.scale_range_low = 0.5
    
    if "ScaleRangeUp" in param:
      self.scale_range_up = float(param["ScaleRangeUp"])
    else:
      self.scale_range_up = 1.5
    
    # ===================Shift=============================
    if "ShiftRangeLow" in param:
      self.shift_range_low = float(param["ShiftRangeLow"])
    else:
      self.shift_range_low = -1.0
    
    if "ShiftRangeUp" in param:
      self.shift_range_up = float(param["ShiftRangeUp"])
    else:
      self.shift_range_up = 1.0
    
    # ===================Syserr=============================
    if "SyserrRangeLow" in param:
      self.syserr_range_low = float(param["SyserrRangeLow"])
    else:
      self.syserr_range_low = 0.0
    
    if "SyserrRangeUp" in param:
      self.syserr_range_up = float(param["SyserrRangeUp"])
    else:
      self.syserr_range_up = 0.1
    
    # ===================Errscale=============================
    if "ErrscaleRangeLow" in param:
      self.errscale_range_low = float(param["ErrscaleRangeLow"])
    else:
      self.errscale_range_low = 0.1
    
    if "ErrscaleRangeUp" in param:
      self.errscale_range_up = float(param["ErrscaleRangeUp"])
    else:
      self.errscale_range_up = 2.0
    

def load_pycali_data_flag(fname):
  """
  load data from a formatted file
  return a dict
  """
  data = {}
  flag = {}
  num_flag = {}
  block = []
  block_flag = []
  flag_list = []
  nc = 0
  code = ""
  fp = open(fname)

  for line in fp:
    if line[0] == "#":
      code = line[1:].split()[0]
      nc = int(line[1:].split()[1])
      block.clear()
      flag_list.clear()
      block_flag.clear()
    else:
      lsp = line.split()

      if len(lsp) == 3:
        block.append(np.array(line.split(), dtype=float))
        flg = 0
      else:
        block.append(np.array(lsp[:3], dtype=float))
        flg = int(lsp[3])

      # flag
      if flg in flag_list:
        idx = np.where(np.array(flag_list) == flg)[0]
        block_flag.append(idx[0])
      else:
        block_flag.append(len(flag_list))
        flag_list.append(flg)

    if len(block) == nc:
      data[code] = np.array(block)
      flag[code] = np.array(block_flag)
      num_flag[code] = len(flag_list)

  fp.close()
  return data, flag, num_flag

def simple_plot(cfg):
  """
  a simple plot
  """ 
  data={}
  nax = 1
  cont = np.loadtxt(cfg.fcont)
  cont_cali = np.loadtxt(cfg.fcont+"_cali", usecols=(0, 1, 2))
  data["cont"]=[cont, cont_cali]
  
  for i, fl in enumerate(cfg.fline):
    nax+=1
    line = np.loadtxt(fl)
    line_cali = np.loadtxt(fl+"_cali", usecols=(0, 1, 2))
    data["line%d"%i] = [line, line_cali]
  
  fig = plt.figure(figsize=(10, 8))
  cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
  for i, key in enumerate(data.keys()):
    ax = fig.add_subplot(nax, 1, i+1)
    d = data[key][0]
    dc = data[key][1]
    ax.errorbar(d[:, 0], d[:, 1], yerr=d[:, 2], ls='none', marker='o', markersize=4, color=cycle[0], 
                ecolor='darkgrey', markeredgecolor=None, elinewidth=1, label=key, fillstyle='none', capsize=0.9)
    ax.errorbar(dc[:, 0], dc[:, 1], yerr=dc[:, 2], ls='none', marker='o', markersize=4, color=cycle[1],
                ecolor='darkgrey', markeredgecolor=None,  elinewidth=1, label=key+r"\_cali", capsize=0.9)
    ax.legend()
    ax.set_xlabel("Time")
    ax.set_ylabel("Flux")
    
    ax.minorticks_on()
  
  plt.show()

def plot_results(cfg):
  """
  a detailed plot, output results to PyCALI_results.pdf.
  """
  pdf = PdfPages("PyCALI_results.pdf")
  
  # whether smooth the histograms
  smooth2d = True
  smooth1d = False

  file_dir = "data"
  file_dir += "/"
  #===================================================================
  # load params
  #===================================================================
  with open(file_dir + "/param_input") as f:
    file_content = '[dump]\n' + f.read()
    
  config = cfgpars.ConfigParser(delimiters='=', allow_no_value=True)
  config.read_string(file_content)
  
  print("================================")
  for key in config["dump"].keys():
    print("%s = %s"%(key, config["dump"][key]))
  print("================================")

  #===================================================================
  # load codes
  #===================================================================
  data, flag, num_flag_cont = load_pycali_data_flag(cfg.fcont)
  # code = np.genfromtxt(file_dir + "/factor.txt", usecols=(0), skip_header=3, dtype=str)
  code = list(data.keys())
  ncode = len(code)
  
  # number of flags
  nsyserr_flag = 0
  for key in data.keys():
    if nsyserr_flag < num_flag_cont[key]:
      nsyserr_flag = num_flag_cont[key]
  
  num_flag_line = []
  for j in range(len(cfg.fline)):
    data, flag, num_flag = load_pycali_data_flag(cfg.fline[j])
    num_flag_line.append(num_flag)
    for key in data.keys():
      if nsyserr_flag < num_flag[key]:
        nsyserr_flag = num_flag[key]

  # remove "_" in code, used for plotting labels
  code_tex = []
  for i in range(len(code)):
    code_tex.append(code[i].replace("_", ""))
  print(code, code_tex)
  
  #===================================================================
  # load means
  #===================================================================
  fp = open(file_dir + "/PyCALI_output.txt", "r")
  cont_mean = np.zeros(ncode)
  fp.readline()
  for i in range(ncode):
    line = fp.readline()
    cont_mean[i] = float(line.split()[2])
  
  lines_mean = {}
  for j in range(len(cfg.fline)):
    lines_mean["%d"%j] = np.zeros(ncode)
    fp.readline()
    for i in range(ncode):
      line = fp.readline()
      lines_mean["%d"%j][i] = float(line.split()[2])
  
  fp.close()
  
  #===================================================================
  # obtain norm for cont and line
  #===================================================================
  num_params_var = 2
  nset = 1
  for j in range(len(cfg.fline)):
    num_params_var += 2
    nset += 1

  cont_mean_code = np.zeros(ncode)
  lines_mean_code={}
  for j in range(len(cfg.fline)):
    line_mean_code = np.zeros(ncode)
    lines_mean_code["%d"%j] = lines_mean_code
  
  sample = np.loadtxt(file_dir + "/posterior_sample.txt")
  # take into account continuum normalization
  sample[:, 0] += np.log(cont_mean[0]) 
  sample[:, 0] /= np.log(10.0) # sigma
  sample[:, 1] /= np.log(10.0) # tau
  for i in range(ncode):
    # scale
    sample[:, num_params_var+i] *= cont_mean[0]/cont_mean[i]
  
    # syserr
    for k in range(nsyserr_flag):
      sample[:, num_params_var+2*ncode + i*nsyserr_flag + k] *= cont_mean[i] 
  
    # error scale
    for k in range(nsyserr_flag):
      sample[:, num_params_var+2*ncode + ncode*nsyserr_flag + i*nsyserr_flag + k] *= 1.0
  
  # shift
  sample[:, num_params_var+ncode:num_params_var+2*ncode] *= cont_mean[0] 
  
  # take into account line normalization
  for j in range(len(cfg.fline)):
    line_mean = lines_mean["%d"%j]
    sample[:, 2+2*j] += np.log(line_mean[0])
    sample[:, 2+2*j] /= np.log(10.0)  # sigma
    sample[:, 2+2*j+1] /= np.log(10.0) # tau

    for i in range(ncode):
      # syserr 
      for k in range(nsyserr_flag):
        sample[:, num_params_var+2*ncode + (2+2*j)*ncode*nsyserr_flag + i*nsyserr_flag + k] *= line_mean[i]
  
      # error scale
      for k in range(nsyserr_flag):
        sample[:, num_params_var+2*ncode + (3+2*j)*ncode*nsyserr_flag + i*nsyserr_flag + k] *=  1.0
  
  # scale in log10
  sample[:, num_params_var:num_params_var+ncode] = np.log10( sample[:, num_params_var:num_params_var+ncode] )
  # error scale in log10
  idx_par = num_params_var+2*ncode+ncode*nsyserr_flag
  sample[:, idx_par:idx_par+ncode*nsyserr_flag] = np.log10( sample[:, idx_par:idx_par+ncode*nsyserr_flag] )
  
  for j in range(len(cfg.fline)):
    idx_par = num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag
    sample[:, idx_par:idx_par+ncode*nsyserr_flag] = np.log10( sample[:, idx_par:idx_par+ncode*nsyserr_flag] )
  
  #===================================================================
  # print posterior values
  #===================================================================
  print("68.3% posterior confidence intervals")
  print("# Note normalization of light curves included,\n# therefore different with those in factor.txt.")
  print("log10 Scale")
  scale = np.zeros(ncode)
  for i in range(ncode):
    mean, low, up = np.quantile(sample[:, num_params_var+i], q=(0.5, 0.16, 0.84))
    scale[i] = mean
    print(code[i], "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
  
  print("\nShift")
  for i in range(ncode):
    mean, low, up = np.quantile(sample[:, num_params_var+ncode+i], q=(0.5, 0.16, 0.84))
    print(code[i], "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
    
  print("\nSyserr of continuum")
  for i in range(ncode):
    for k in range(nsyserr_flag):
      mean, low, up = np.quantile(sample[:, num_params_var+2*ncode+i*nsyserr_flag+k], q=(0.5, 0.16, 0.84))
      print(code[i], k, "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
  
  print("\nlog10 Error Scale of continuum")
  for i in range(ncode):
    for k in range(nsyserr_flag):
      mean, low, up = np.quantile(sample[:, num_params_var+2*ncode+ncode*nsyserr_flag+i*nsyserr_flag+k], q=(0.5, 0.16, 0.84))
      print(code[i], k, "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
  
  for j in range(len(cfg.fline)):
    print("\nSyserr of line%d"%j)
    for i in range(ncode):
      for k in range(nsyserr_flag):
        mean, low, up = np.quantile(sample[:, num_params_var+2*ncode+(2+2*j)*ncode*nsyserr_flag+i*nsyserr_flag+k], q=(0.5, 0.16, 0.84))
        print(code[i], k, "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
  
    print("\nlog10 Error Scale of line%d"%j)
    for i in range(ncode):
      for k in range(nsyserr_flag):
        mean, low, up = np.quantile(sample[:, num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag+i*nsyserr_flag+k], q=(0.5, 0.16, 0.84))
        print(code[i], k, "%5.3e -%5.3e +%5.3e"%(mean, mean-low, up-mean))
  print("================================")
  #exit()

  plt.rc('text', usetex=True)
  plt.rc('font', family="serif", size=18)
  
  #===================================================================
  # now plot
  #===================================================================
  print("ploting, wait...")
  data={}
  nax = 1
  # first continuum 
  cont = np.loadtxt(cfg.fcont)
  cont_code_org = np.empty(cont.shape[0], dtype="U20")
  cont_cali = np.loadtxt(cfg.fcont+"_cali", usecols=(0, 1, 2))
  cont_code = np.loadtxt(cfg.fcont+"_cali", usecols=(3), dtype=str)
  data["cont"]=[cont, cont_cali]
  cont_full = np.loadtxt(cfg.fcont+"_recon")
  
  # create original code of the raw data
  i1=0
  i2=0
  for i in range(ncode):
    nc = np.count_nonzero(cont_code==code[i])
    i2 = i1 + nc
    cont_code_org[i1:i2]=code[i]

    if nc > 0:
      cont_mean_code[i] = np.mean(cont[i1:i2, 1])
    else:
      cont_mean_code[i] = 1.0
    i1 = i2
  
  # load index for sorting the data
  idx_cont = np.loadtxt(cfg.fcont+"_sort", dtype=int)
  
  # load line data if included
  lines_code = {}
  lines_code_org = {}
  lines_full = {}
  idx_lines = {}
  for j in range(len(cfg.fline)):
    line = np.loadtxt(cfg.fline[j])
    line_code_org = np.empty(line.shape[0], dtype="U20")
    line_cali = np.loadtxt(cfg.fline[j]+"_cali", usecols=(0, 1, 2))
    line_code = np.loadtxt(cfg.fline[j]+"_cali", usecols=(3), dtype=str)
    data["line%d"%j] = [line, line_cali]
    line_full = np.loadtxt(cfg.fline[j]+"_recon")
    idx_line = np.loadtxt(cfg.fline[j]+"_sort", dtype=int)
    
    i1=0
    i2=0
    line_mean_code = lines_mean_code["%d"%j]
    for i in range(ncode):
      nc = np.count_nonzero(line_code==code[i])
      i2 = i1 + nc
      line_code_org[i1:i2]=code[i]

      if nc > 0:
        line_mean_code[i] = np.mean(line[i1:i2, 1])
      else:
        line_mean_code[i] = 1.0

      i1 = i2
    
    lines_code["%d"%j] = line_code
    lines_code_org["%d"%j] = line_code_org
    lines_full["%d"%j] = line_full
    idx_lines["%d"%j] = idx_line
  
  # obtain colors of matplotlib
  if ncode <= 10: 
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
  else:
    cycle = [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
        '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
        '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
        '#17becf', '#9edae5']
  
  fig = plt.figure(figsize=(15, 12))
  
  # plot original data
  ax = fig.add_axes((0.1, 0.68, 0.66, 0.28))
  key="cont"
  d = data[key][0]
  dc = data[key][1]
  for i in range(ncode):
    idx = np.where((cont_code_org == code[i]))
    ax.errorbar(d[idx[0], 0], d[idx[0], 1], yerr=d[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                ecolor='grey', markeredgecolor=None, elinewidth=1, capsize=0.9,  label=r'{0} $\rm {1}~({2})$'.format(i, code_tex[i], len(idx[0])))
                
  ax.legend(frameon=False, loc=(1.0, 0.0), handletextpad=-0.1, fontsize=15)
  ax.set_ylabel("Raw Data Flux")
  xlim = ax.get_xlim()
  ylim = ax.get_ylim()
  [xt.set_visible(False) for xt in ax.get_xticklabels()]
  ax.minorticks_on()
  
  # plot calibrated data              
  ax = fig.add_axes((0.1, 0.38, 0.66, 0.28))
  ax.plot(cont_full[:, 0], cont_full[:, 1], lw=1, linestyle="--", color='k', alpha=0.8)
  ax.plot(cont_full[:, 0], cont_full[:, 1]+cont_full[:, 2], lw=1, linestyle="--", color='k', alpha=0.8)
  ax.plot(cont_full[:, 0], cont_full[:, 1]-cont_full[:, 2], lw=1, linestyle="--", color='k', alpha=0.8)
  for i in range(ncode):
    idx = np.where((cont_code == code[i]))
    # errors after calibration
    ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=1.5)
    # errors only from scaling
    ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=d[idx_cont[idx[0]], 2]*10**(scale[i]), ls='none', color=cycle[np.mod(i, len(cycle), dtype=int)], \
                ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=1.5)
    
  ax.set_ylabel("Intercalibrated Flux")
  ax.set_xlim(xlim[0], xlim[1])
  [xt.set_visible(False) for xt in ax.get_xticklabels()]
  ax.minorticks_on()
  
  # plot parameter prior
  ax = fig.add_axes((0.76, 0.38, 0.2, 0.5))
  ax.text(0.3, 0.5, r"$\varphi,~~~~G, ~~~\epsilon, ~~~b$", fontsize=15)
  for i in range(ncode):
    fstr = r"${0}$".format(i)
    ax.text(0.1, 0.45-i*0.04, fstr, fontsize=15)
    fstr=r""
    if np.std(sample[:, num_params_var+i]) <= EPS :
      fstr = fstr + r"N"
    else:
      fstr = fstr + r"Y"
    
    if np.std(sample[:, num_params_var+ncode+i]) <= EPS :
      fstr = fstr + r"~~~~~N"
    else:
      fstr = fstr + r"~~~~~Y"
    
    if np.std(sample[:, num_params_var+2*ncode+i*nsyserr_flag]) <= EPS :
      fstr = fstr + r"~~~~N"
    else:
      fstr = fstr + r"~~~~Y"
    
    if np.std(sample[:, num_params_var+2*ncode+ncode*nsyserr_flag+i*nsyserr_flag]) <= EPS :
      fstr = fstr + r"~~~N"
    else:
      fstr = fstr + r"~~~Y"
    
    ax.text(0.3, 0.45-i*0.04, fstr, fontsize=15)
  
  ax.text(0.1, 0.45-ncode*0.04, "Y: free, N: fixed", fontsize=15)
  ax.set_axis_off()
  
  # plot residuals
  ax = fig.add_axes((0.1, 0.08, 0.66, 0.28))
  for i in range(ncode):
   idx = np.where((cont_code == code[i]))
   res = dc[idx[0], 1] - np.interp(dc[idx[0], 0], cont_full[:, 0], cont_full[:, 1])
   # errors after calibration
   ax.errorbar(dc[idx[0], 0], res, yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                ecolor=cycle[np.mod(i, len(cycle))], markeredgecolor=None,  elinewidth=1, capsize=1.5, zorder=1)
   # errors only from scaling
   ax.errorbar(dc[idx[0], 0], res, yerr=d[idx_cont[idx[0]], 2]*10**(scale[i]), ls='none', color=cycle[np.mod(i, len(cycle))], \
                ecolor=cycle[np.mod(i, len(cycle))], markeredgecolor=None,  elinewidth=1, capsize=1.5, zorder=1)
  
  error_mean = np.mean(dc[:, 2])
  ax.axhline(y=0.0, linestyle='--', color='silver', lw=1, zorder=0)
  ax.axhline(y=-error_mean, linestyle='--', color='silver', lw=1, zorder=0)
  ax.axhline(y=error_mean, linestyle='--', color='silver', lw=1, zorder=0)
  ax.set_xlabel("Time")
  #ax.set_title("Continuum Residuals")
  ax.set_ylabel("Residuals")
  ax.set_xlim(xlim[0], xlim[1])
  ylim = ax.get_ylim()
  ax.minorticks_on()
  
  ax = fig.add_axes((0.76, 0.08, 0.07, 0.28))
  xlim = ax.get_xlim()
  for i in range(ncode):
    idx = np.where((cont_code == code[i]))
    if len(idx[0]) == 0:
      continue
    # errors after calibration
    ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(dc[idx[0], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
               elinewidth=1, capsize=1.5, zorder=1)
    # errors only from scaling
    ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(d[idx_cont[idx[0]], 2])*10**(scale[i]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
               elinewidth=1, capsize=1.5, zorder=1)
  
  ax.set_xlim(xlim[0], xlim[1])
  ax.set_ylim(ylim[0], ylim[1])
  [xt.set_visible(False) for xt in ax.get_xticklabels()]
  [xt.set_visible(False) for xt in ax.get_yticklabels()]
  ax.minorticks_on()
  ax.axhline(y=0.0, linestyle='--', color='silver', lw=1, zorder=0)
  ax.axhline(y=-error_mean, linestyle='--', color='silver', lw=1, zorder=0)
  ax.axhline(y=error_mean, linestyle='--', color='silver', lw=1, zorder=0)
  
  ax = fig.add_axes((0.88, 0.08, 0.1, 0.28))
  ax.hist((dc[:, 1] - np.interp(dc[:, 0], cont_full[:, 0], cont_full[:, 1]))/dc[:, 2], orientation='horizontal', \
         density=True, bins=20, range=[-4, 4])
  y = np.linspace(-4, 4, 100)
  x = 1.0/np.sqrt(2.0*np.pi)*np.exp(-0.5*y*y)
  ax.plot(x, y)
  ax.set_ylim(-4, 4)
  #[yt.set_visible(False) for yt in ax.get_yticklabels()]
  ax.set_ylabel("Stardarized Residuals")
  ax.minorticks_on()
  
  fname = cfg.fcont
  fname = fname.replace("_", "\_")
  fig.suptitle(r"\bf {0}".format(fname), x=0.5, y=1.0)
  pdf.savefig(fig)
  plt.close()
  
  #===================================================================
  # then plot line if there is
  #===================================================================
  for j in range(len(cfg.fline)):
    fig = plt.figure(figsize=(15, 12))
    
    ax = fig.add_axes((0.1, 0.68, 0.66, 0.28))
    key="line%d"%j
    d = data[key][0]
    dc = data[key][1]
    line_code_org = lines_code_org["%d"%j]
    line_code = lines_code["%d"%j]
    line_full = lines_full["%d"%j]
    idx_line = idx_lines["%d"%j]
    for i in range(ncode):
     idx = np.where((line_code_org == code[i]))
     ax.errorbar(d[idx[0], 0], d[idx[0], 1], yerr=d[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                 ecolor='grey', markeredgecolor=None, elinewidth=1, capsize=0.9,  label=r'{0} $\rm {1}~({2})$'.format(i, code_tex[i], len(idx[0])))
    
    ax.legend(frameon=False, loc=(1.0, 0.0), handletextpad=-0.1, fontsize=15)
    ax.set_ylabel("Raw Data Flux")
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.minorticks_on()
    [xt.set_visible(False) for xt in ax.get_xticklabels()]
   
    ax = fig.add_axes((0.1, 0.38, 0.66, 0.28))
    ax.plot(line_full[:, 0], line_full[:, 1], lw=1, linestyle="--", color='k', alpha=0.8)
    ax.plot(line_full[:, 0], line_full[:, 1]-line_full[:, 2], lw=1, linestyle="--", color='k', alpha=0.8)
    ax.plot(line_full[:, 0], line_full[:, 1]+line_full[:, 2], lw=1, linestyle="--", color='k', alpha=0.8)
    for i in range(ncode):
      idx = np.where((line_code == code[i]))
      # errors after calibration
      ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                 ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=0.9, label=r'${0}$'.format(code_tex[i]))
      
      # errors only from scaling
      ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=d[idx_line[idx[0]], 2]*10**(scale[i]), ls='none', color=cycle[np.mod(i, len(cycle), dtype=int)], \
                ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=1.5)
    
    ax.set_ylabel("Intercalibrated Flux")
    ax.set_xlim(xlim[0], xlim[1])
    ax.minorticks_on()
    [xt.set_visible(False) for xt in ax.get_xticklabels()]
   
    # plot parameter prior
    ax = fig.add_axes((0.76, 0.38, 0.2, 0.5))
    ax.text(0.3, 0.5, r"$\varphi,~~~~G, ~~~\epsilon, ~~~b$", fontsize=15)
    for i in range(ncode):
      fstr = r"${0}$".format(i)
      ax.text(0.1, 0.45-i*0.04, fstr, fontsize=15)
      fstr=r""
      if np.std(sample[:, num_params_var+i]) <= EPS :
        fstr = fstr + r"N"
      else:
        fstr = fstr + r"Y"
      
      # line does not have G
      if np.std(sample[:, num_params_var+ncode+i]) <= EPS :
        fstr = fstr + r"~~~~~N"
      else:
        fstr = fstr + r"~~~~~N"
      
      if np.std(sample[:, num_params_var+2*ncode+i*nsyserr_flag]) <= EPS :
        fstr = fstr + r"~~~~N"
      else:
        fstr = fstr + r"~~~~Y"
      
      if np.std(sample[:, num_params_var+2*ncode+ncode*nsyserr_flag+i*nsyserr_flag]) <= EPS :
        fstr = fstr + r"~~~N"
      else:
        fstr = fstr + r"~~~Y"
      
      ax.text(0.3, 0.45-i*0.04, fstr, fontsize=15)
   
    ax.text(0.1, 0.45-ncode*0.04, "Y: free, N: fixed", fontsize=15)
    ax.set_axis_off()
   
    ax = fig.add_axes((0.1, 0.08, 0.66, 0.28))
    for i in range(ncode):
      idx = np.where((line_code == code[i]))
      res = dc[idx[0], 1] - np.interp(dc[idx[0], 0], line_full[:, 0], line_full[:, 1])
      # errors after calibration
      ax.errorbar(dc[idx[0], 0], res, yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
                 ecolor=cycle[np.mod(i, len(cycle), dtype=int)], markeredgecolor=None,  elinewidth=1, capsize=0.9, label=r'${0}$'.format(code_tex[i]), zorder=0)
      
      # errors only from scaling 
      ax.errorbar(dc[idx[0], 0], res, yerr=d[idx_line[idx[0]], 2]*10**(scale[i]), \
                  ls='none', color=cycle[np.mod(i, len(cycle))], \
                  ecolor=cycle[np.mod(i, len(cycle))], markeredgecolor=None,  elinewidth=1, capsize=1.5, zorder=1)
      
    error_mean = np.mean(dc[:, 2])
    ax.axhline(y=0.0, linestyle='--', color='silver', lw=1, zorder=0)
    ax.axhline(y=-error_mean, linestyle='--', color='silver', lw=1, zorder=0)
    ax.axhline(y=error_mean, linestyle='--', color='silver', lw=1, zorder=0)
    ax.set_xlabel("Time")
    ax.set_ylabel("Residuals")
    ax.set_xlim(xlim[0], xlim[1])
    ylim = ax.get_ylim()
    ax.minorticks_on()
    
    ax = fig.add_axes((0.76, 0.08, 0.07, 0.28))
    xlim = ax.get_xlim()
    for i in range(ncode):
      idx = np.where((line_code == code[i]))
      if len(idx[0]) == 0:
        continue
      # errors after calibration
      ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(dc[idx[0], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
                 elinewidth=1, capsize=1.5, zorder=1)
      # errors only from scaling
      ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(d[idx_line[idx[0]], 2])*10**(scale[i]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
                 elinewidth=1, capsize=1.5, zorder=1)
    
    ax.set_xlim(xlim[0], xlim[1])
    ax.set_ylim(ylim[0], ylim[1])
    [xt.set_visible(False) for xt in ax.get_xticklabels()]
    [xt.set_visible(False) for xt in ax.get_yticklabels()]
    ax.minorticks_on()
    ax.axhline(y=0.0, linestyle='--', color='silver', lw=1, zorder=0)
    ax.axhline(y=-error_mean, linestyle='--', color='silver', lw=1, zorder=0)
    ax.axhline(y=error_mean, linestyle='--', color='silver', lw=1, zorder=0)
   
    ax = fig.add_axes((0.88, 0.08, 0.1, 0.28))
    ax.hist((dc[:, 1] - np.interp(dc[:, 0], line_full[:, 0], line_full[:, 1]))/dc[:, 2], orientation='horizontal', \
            density=True, bins=20, range=[-4, 4])
    y = np.linspace(-4, 4, 100)
    x = 1.0/np.sqrt(2.0*np.pi)*np.exp(-0.5*y*y)
    ax.plot(x, y)
    ax.set_ylim(-4, 4)
    
    #[yt.set_visible(False) for yt in ax.get_yticklabels()]
    ax.set_ylabel("Stardarized Residuals")
    ax.minorticks_on()
    
    fname = cfg.fline[j]
    fname = fname.replace("_", "\_")
    fig.suptitle(r"\bf {0}".format(fname), x=0.5, y=1.0)
   
    pdf.savefig(fig)
    plt.close()
   
  #pdf.close()
  #exit()
  
  #===================================================================
  # now plot histograms
  #===================================================================
  
  # first sigma and tau
  fig = corner.corner(sample[:, :2], smooth=smooth2d, smooth1d = smooth1d, labels=[r"$\log\sigma$", r"$\log\tau$"], \
        levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
  
  axes = np.array(fig.axes).reshape((2, 2))
  ax = axes[0, 1]
  ax.text(0.0, 0.5, r"\bf $\sigma$ \& $\tau$", fontsize=15)
  ax.text(0.0, 0.65, r"\bf cont", fontsize=15)
  # plot limits
  # sigma
  ax = axes[0, 0]
  xlim = ax.get_xlim()
  if(xlim[0]<np.log10(cfg.sigma_range_low * cont_mean[0])):
    ax.axvline(x=np.log10(cfg.sigma_range_low * cont_mean[0]), ls='--')
  if xlim[1]>np.log10(cfg.sigma_range_up * cont_mean[0]):
    ax.axvline(x=np.log10(cfg.sigma_range_up * cont_mean[0]), ls='--')
  
  # tau
  ax = axes[1, 1]
  xlim = ax.get_xlim()
  if(xlim[0]<np.log10(cfg.tau_range_low)):
    ax.axvline(x=np.log10(cfg.tau_range_low), ls='--')
  if xlim[1]>np.log10(cfg.tau_range_up):
    ax.axvline(x=np.log10(cfg.tau_range_up), ls='--')
  pdf.savefig(fig)
  plt.close()

  for j in range(len(cfg.fline)):
    line_mean = lines_mean["%d"%j]
    fig = corner.corner(sample[:, (j+1)*2:(j+2)*2], smooth=smooth2d, smooth1d = smooth1d, labels=[r"$\log\sigma$", r"$\log\tau$"], \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
  
    axes = np.array(fig.axes).reshape((2, 2))
    ax = axes[0, 1]
    ax.text(0.0, 0.5, r"\bf $\sigma$ \& $\tau$", fontsize=15)
    ax.text(0.0, 0.65, r"\bf"+r'$\bf line {0}$'.format(j), fontsize=15)

    # plot limits
    # sigma
    ax = axes[0, 0]
    xlim = ax.get_xlim()
    if(xlim[0]<np.log10(cfg.sigma_range_low * line_mean[i])):
      ax.axvline(x=np.log10(cfg.sigma_range_low * line_mean[i]), ls='--')
    if(xlim[1]>np.log10(cfg.sigma_range_up * line_mean[i])):
      ax.axvline(x=np.log10(cfg.sigma_range_up * line_mean[i]), ls='--')
    
    # tau
    ax = axes[1, 1]
    xlim = ax.get_xlim()
    if(xlim[0]<np.log10(cfg.tau_range_low)):
      ax.axvline(x=np.log10(cfg.tau_range_low), ls='--')
    if xlim[1]>np.log10(cfg.tau_range_up):
      ax.axvline(x=np.log10(cfg.tau_range_up), ls='--')
    
    pdf.savefig(fig)
    plt.close()


  if int(config["dump"]["fixed_scale"]) == 1:
    fig = corner.corner(sample[:, num_params_var+ncode+1:num_params_var+2*ncode], smooth=smooth2d, smooth1d = smooth1d, \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f", range=[0.9999]*(ncode-1))

    axes = fig.get_axes()
    for i in range(ncode-1):
      ax = axes[i*(ncode-1)+i]
      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i+1]))

      # plot limits
      if(xlim[0]<cfg.shift_range_low * cont_mean[0]):
        ax.axvline(x=cfg.shift_range_low * cont_mean[0], ls='--')
      if xlim[1]>cfg.shift_range_up * cont_mean[0]:
        ax.axvline(x=cfg.shift_range_up * cont_mean[0], ls='--')
      

    fig.suptitle(r"\bf Shift", fontsize=20)
    pdf.savefig(fig)
    plt.close()
   
  elif int(config["dump"]["fixed_shift"]) == 1:
    fig = corner.corner(sample[:, num_params_var+1:num_params_var+ncode], smooth=smooth2d, smooth1d = smooth1d,  \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f", range=[0.9999]*(ncode-1))
    axes = fig.get_axes()
    for i in range(ncode-1):
      ax = axes[i*(ncode-1)+i]
      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i+1]))
      
      # plot limits
      if(xlim[0]<np.log10(cfg.scale_range_low * cont_mean[0]/cont_mean[i+1])):
        ax.axvline(x=np.log10(cfg.scale_range_low * cont_mean[0]/cont_mean[i+1]), ls='--')
      if xlim[1]>np.log10(cfg.scale_range_up * cont_mean[0]/cont_mean[i+1]):
        ax.axvline(x=np.log10(cfg.scale_range_up * cont_mean[0]/cont_mean[i+1]), ls='--')
      

    fig.suptitle(r"\bf Scale", fontsize=20)
    pdf.savefig(fig)
    plt.close()
  else:
    for i in range(1, ncode):
      range_min = np.min(sample[:, [num_params_var+i,num_params_var+i+ncode]], axis=0)
      range_max = np.max(sample[:, [num_params_var+i,num_params_var+i+ncode]], axis=0)
      span = range_max - range_min
      range_interval = [[range_min[i]-0.3*span[i], range_max[i]+0.3*span[i]] for i in range(2)]
      fig = corner.corner(sample[:, [num_params_var+i,num_params_var+i+ncode]], smooth=smooth2d, smooth1d = smooth1d, labels=[r"$\log\varphi$", r"$G$"], 
            range=range_interval, \
            levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
      
      # plot input limits
      axes = np.array(fig.axes).reshape((2, 2))
      # scale
      ax = axes[0, 0]
      xlim = ax.get_xlim()
      if(xlim[0]<np.log10(cfg.scale_range_low * cont_mean[0]/cont_mean[i])):
        ax.axvline(x=np.log10(cfg.scale_range_low * cont_mean[0]/cont_mean[i]), ls='--')
      if(xlim[1]>np.log10(cfg.scale_range_up * cont_mean[0]/cont_mean[i])):
        ax.axvline(x=np.log10(cfg.scale_range_up * cont_mean[0]/cont_mean[i]), ls='--')
      
      # shift
      ax = axes[1, 1]
      xlim = ax.get_xlim()
      if(xlim[0]<cfg.shift_range_low * cont_mean[0]):
        ax.axvline(x=cfg.shift_range_low * cont_mean[0], ls='--')
      if xlim[1]>cfg.shift_range_up * cont_mean[0]:
        ax.axvline(x=cfg.shift_range_up * cont_mean[0], ls='--')

      ax = axes[0, 1]
      ax.text(0.0, 0.5, r"\bf Scale \& Shift", fontsize=15)
      ax.text(0.0, 0.65, r"\bf"+r'$\bf {0}$'.format(code_tex[i]), fontsize=15)
      #fig.suptitle(r"\bf Scale \& Shift  "+r'${0}$'.format(code[i][3:-4]), fontsize=20)
      pdf.savefig(fig)
      plt.close()
  
  if int(config["dump"]["fixed_syserr"]) == 0:
    if nsyserr_flag == 1:
      idx_par = num_params_var+2*ncode
      fig = corner.corner(sample[:, idx_par:idx_par+ncode*nsyserr_flag], smooth=smooth2d, smooth1d = smooth1d, \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f", range=[0.9999]*ncode)
      axes = fig.get_axes()
      for i in range(ncode):
        ax = axes[i*ncode+i]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i]))

        # plot limits 
        if xlim[0] < cfg.syserr_range_low * cont_mean[i]:
          ax.axvline(x=cfg.syserr_range_low * cont_mean[i], ls='--')
        if xlim[1]>cfg.syserr_range_up * cont_mean[i]:
          ax.axvline(x=cfg.syserr_range_up * cont_mean[i], ls='--')

      fig.suptitle(r"\bf Systematic Error (Cont)", fontsize=20)
      pdf.savefig(fig)
      plt.close()
    else:
      keys = list(num_flag_cont.keys())
      for i in range(ncode):
        if num_flag_cont[keys[i]] == 0:
          continue
        idx_par = num_params_var+2*ncode+i*nsyserr_flag
        fig = corner.corner(sample[:, idx_par:idx_par+num_flag_cont[keys[i]]],
                            smooth=smooth2d, smooth1d = smooth1d, \
                            levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), \
                            show_titles=False, title_fmt=".3f")
        fig.suptitle(r"\bf Systematic Error (Cont %s)"%keys[i], fontsize=12)
        pdf.savefig(fig)
        plt.close() 

  if int(config["dump"]["fixed_error_scale"]) == 0:
    if nsyserr_flag == 1:
      idx_par = num_params_var+2*ncode+ncode
      fig = corner.corner(sample[:, idx_par:idx_par+ncode*nsyserr_flag], smooth=smooth2d, smooth1d = smooth1d, \
          show_titles=True, title_fmt=".3f", range=[0.9999]*ncode)
      axes = fig.get_axes()
      for i in range(ncode):
        ax = axes[i*ncode+i]
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i]))

        # plot limits 
        if xlim[0] < np.log10(cfg.errscale_range_low):
          ax.axvline(x=np.log10(cfg.errscale_range_low), ls='--')
        if xlim[1]>np.log10(cfg.errscale_range_up):
          ax.axvline(x=np.log10(cfg.errscale_range_up), ls='--')

      fig.suptitle(r"\bf log(Error Scale) (Cont)", fontsize=20)
      pdf.savefig(fig)
      plt.close()

    else:
      keys = list(num_flag_cont.keys())
      for i in range(ncode):
        if num_flag_cont[keys[i]] == 0:
          continue

        idx_par = num_params_var+2*ncode+ncode*nsyserr_flag+i*nsyserr_flag
        fig = corner.corner(sample[:, idx_par:idx_par+num_flag_cont[keys[i]]],
                            smooth=smooth2d, smooth1d = smooth1d, \
                            levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), \
                            show_titles=False, title_fmt=".3f")
        fig.suptitle(r"\bf log(Error Scale) (Cont %s)"%keys[i], fontsize=12)
        pdf.savefig(fig)
        plt.close() 
    
  if int(config["dump"]["fixed_syserr"]) == 0 and int(config["dump"]["fixed_error_scale"]) == 0:
    keys = list(num_flag_cont.keys())
    for i in range(ncode):
      for j in range(num_flag_cont[keys[i]]):
        idx_par = num_params_var+2*ncode+i*nsyserr_flag+j
        fig = corner.corner(sample[:, [idx_par,idx_par+ncode*nsyserr_flag]], smooth=smooth2d, smooth1d = smooth1d, labels=[r"$\epsilon$", r"$b$"], 
              levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
        
        axes = fig.get_axes()
        ax = axes[1]
        ax.text(0.0, 0.5, r"\bf Syserr \& log(Error Scale)", fontsize=12)
        ax.text(0.0, 0.6, r"\bf "+r"{0} {1}".format(code_tex[i], j), fontsize=15)
        #fig.suptitle(r"\bf Syserr \& Error Scale  "+code[i][3:-4], fontsize=20)

        # plot limits
        ax = axes[0]
        xlim = ax.get_xlim()
        if xlim[0] < cfg.syserr_range_low * cont_mean[i]:
          ax.axvline(x=cfg.syserr_range_low * cont_mean[i], ls='--')
        if(xlim[1]>cfg.syserr_range_up * cont_mean[i]):
          ax.axvline(x=cfg.syserr_range_up * cont_mean[i], ls='--')
        
        ax = axes[1*2+1]
        xlim = ax.get_xlim()
        if xlim[0] < np.log10(cfg.errscale_range_low):
          ax.axvline(x=np.log10(cfg.errscale_range_low), ls='--')
        if xlim[1]>np.log10(cfg.errscale_range_up):
          ax.axvline(x=np.log10(cfg.errscale_range_up), ls='--')

        pdf.savefig(fig)
        plt.close()
  
  # histograms for line
  for j in range(len(cfg.fline)):
    line_mean = lines_mean["%d"%j]
    if int(config["dump"]["fixed_syserr"]) == 0:
      if nsyserr_flag == 1:
        fig = corner.corner(sample[:, num_params_var+4*ncode+2*j*ncode:num_params_var+5*ncode+2*j*ncode], smooth=smooth2d, smooth1d = smooth1d, \
              levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f", range=[0.9999]*ncode)
        axes = fig.get_axes()
        for i in range(ncode):
          ax = axes[i*ncode+i]
          xlim = ax.get_xlim()
          ylim = ax.get_ylim()
          ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i]))
          
          # plot limits
          if xlim[0] < cfg.syserr_range_low * line_mean[i]:
            ax.axvline(x=cfg.syserr_range_low * line_mean[i], ls='--')
          if xlim[1]>cfg.syserr_range_up * line_mean[i]:
            ax.axvline(x=cfg.syserr_range_up * line_mean[i], ls='--')

        fig.suptitle(r"\bf Systematic Error (Line%d)"%j, fontsize=20)
        pdf.savefig(fig)
        plt.close()
      else:
        num_flag = num_flag_line[j]
        keys = list(num_flag.keys())
        for i in range(ncode):
          if num_flag[keys[i]] == 0:
            continue

          idx_par = num_params_var+2*ncode+(2+2*j)*ncode*nsyserr_flag+i*nsyserr_flag
          fig = corner.corner(sample[:, idx_par:idx_par+num_flag[keys[i]]],
                              smooth=smooth2d, smooth1d = smooth1d, \
                              levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), \
                              show_titles=False, title_fmt=".3f")
          fig.suptitle(r"\bf Systematic Error (Line %d %s)"%(j, keys[i]), fontsize=12)
          pdf.savefig(fig)
          plt.close() 

    if int(config["dump"]["fixed_error_scale"]) == 0:
      if nsyserr_flag == 1:
        fig = corner.corner(sample[:, num_params_var+5*ncode+2*j*ncode:num_params_var+6*ncode+2*j*ncode], smooth=smooth2d, smooth1d = smooth1d,\
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f", range=[0.9999]*ncode)
        axes = fig.get_axes()
        for i in range(ncode):
          ax = axes[i*ncode+i]
          xlim = ax.get_xlim()
          ylim = ax.get_ylim()
          ax.text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code_tex[i]))

          # plot limits
          if xlim[0] < np.log10(cfg.errscale_range_low):
            ax.axvline(x=np.log10(cfg.errscale_range_low), ls='--')
          if xlim[1]>np.log10(cfg.errscale_range_up):
            ax.axvline(x=np.log10(cfg.errscale_range_up), ls='--')

        fig.suptitle(r"\bf log(Error Scale) (Line%d)"%j, fontsize=20)
        pdf.savefig(fig)
        plt.close()
      else:
        num_flag = num_flag_line[j]
        keys = list(num_flag.keys())
        for i in range(ncode):
          if num_flag[keys[i]] == 0:
            continue

          idx_par = num_params_var+2*ncode+(3+2*j)*ncode*nsyserr_flag+i*nsyserr_flag
          fig = corner.corner(sample[:, idx_par:idx_par+num_flag[keys[i]]],
                              smooth=smooth2d, smooth1d = smooth1d, \
                              levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), \
                              show_titles=False, title_fmt=".3f")
          fig.suptitle(r"\bf log(Error Scale) (Line %d %s)"%(j, keys[i]), fontsize=12)
          pdf.savefig(fig)
          plt.close() 
  
    if int(config["dump"]["fixed_syserr"]) == 0 and int(config["dump"]["fixed_error_scale"]) == 0:
      num_flag = num_flag_line[j]
      keys = list(num_flag.keys())
      for i in range(ncode):
        for k in range(num_flag[keys[i]]):
          idx_par = num_params_var+2*ncode+(2+2*j)*ncode*nsyserr_flag+i*nsyserr_flag+k
          fig = corner.corner(sample[:, [idx_par,idx_par+ncode*nsyserr_flag]], smooth=smooth2d, smooth1d = smooth1d, \
              labels=[r"$\epsilon$", r"$b$"], levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
        
          axes = fig.get_axes()
          ax = axes[1]
          ax.text(0.0, 0.5, r"\bf Syserr \& log(Error Scale)", fontsize=12)
          ax.text(0.0, 0.65, r"\bf "+r'{0} {1}'.format(code_tex[i], k), fontsize=12)

          # plot limits
          ax = axes[0]
          xlim = ax.get_xlim()
          if xlim[0] < cfg.syserr_range_low * line_mean[i]:
            ax.axvline(x=cfg.syserr_range_low * line_mean[i], ls='--')
          if xlim[1]>cfg.syserr_range_up * line_mean[i]:
            ax.axvline(x=cfg.syserr_range_up * line_mean[i], ls='--')
          
          ax = axes[1*2+1]
          xlim = ax.get_xlim()
          if xlim[0] < np.log10(cfg.errscale_range_low):
            ax.axvline(x=np.log10(cfg.errscale_range_low), ls='--')
          if xlim[1]>np.log10(cfg.errscale_range_up):
            ax.axvline(x=np.log10(cfg.errscale_range_up), ls='--')

          pdf.savefig(fig)
          plt.close()
  
    
  pdf.close()

if __name__=="__main__":
  if len(sys.argv) < 2:
    raise Exception("Please specify a param file as: 'python plot_for_cali.py param.txt'.")

  # load configuration form param.txt
  cfg = Config(sys.argv[1])
   
  # plot results to PyCALI_results.pdf
  plot_results(cfg)
  
  # a simple plot 
  simple_plot(cfg)