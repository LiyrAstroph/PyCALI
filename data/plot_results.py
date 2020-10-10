#
# plot hists for posterior samples generated by pyCali
#
import numpy as np
import matplotlib.pyplot as plt 
import corner
import configparser as cfgpars
from matplotlib.backends.backend_pdf import PdfPages
from os.path import basename
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

pdf = PdfPages("pyCALI_results.pdf")

#===================================================================
# load params
#===================================================================
with open("param_input") as f:
  file_content = '[dump]\n' + f.read()
  
config = cfgpars.ConfigParser(delimiters='=', allow_no_value=True)
config.read_string(file_content)

for key in config["dump"].keys():
  print(key, config["dump"][key])

#===================================================================
# obtain norm for cont and line
#===================================================================
cont = np.loadtxt(basename(config["dump"]["fcont"]))
fp=open(basename(config["dump"]["fcont"]), "r")
fstr=fp.readline()
fstr=fstr[1:].strip()
num = int(fstr.split()[1])
norm_cont = np.mean(cont[:num, 1])
fp.close()

#===================================================================
# load codes
#===================================================================
code = np.genfromtxt("factor.txt", usecols=(0), skip_header=1, dtype=str)

num_params_var = 2
nset = 1
if config["dump"]["fline"] != "":
  num_params_var += 2
  nset += 1
  line = np.loadtxt(basename(config["dump"]["fline"]))
  fp=open(basename(config["dump"]["fline"]), "r")
  fstr=fp.readline()
  fstr=fstr[1:].strip()
  num = int(fstr.split()[1])
  norm_line = np.mean(line[:num, 1])
  fp.close()
 
sample = np.loadtxt("posterior_sample.txt")
ncode = (sample.shape[1] - num_params_var)//(2+nset*2)

# take into account continuum normalization
sample[:, 0] += np.log(norm_cont) 
sample[:, 0] /= np.log(10.0)
sample[:, num_params_var+ncode:num_params_var+2*ncode] *= norm_cont 
sample[:, num_params_var+2*ncode:num_params_var+3*ncode] *= norm_cont

# take into account line normalization
if config["dump"]["fline"] != "":
  sample[:, 2] += np.log(norm_line)
  sample[:, 2] /= np.log(10.0)
  sample[:, num_params_var+4*ncode:num_params_var+5*ncode] *= norm_line

# scale in log10
sample[:, num_params_var:num_params_var+ncode] = np.log10( sample[:, num_params_var:num_params_var+ncode] )
sample[:, num_params_var+3*ncode:num_params_var+4*ncode] = np.log10( sample[:, num_params_var+3*ncode:num_params_var+4*ncode] )
if config["dump"]["fline"] != "":
  sample[:, num_params_var+5*ncode:num_params_var+6*ncode] = np.log10( sample[:, num_params_var+5*ncode:num_params_var+6*ncode] )

#===================================================================
# print posterior values
#===================================================================
print("log10 Scale")
for i in range(ncode):
  mean, low, up = np.quantile(sample[:, num_params_var+i], q=(0.5, 0.16, 0.84))
  print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))

print("\nShift")
for i in range(ncode):
  mean, low, up = np.quantile(sample[:, num_params_var+ncode+i], q=(0.5, 0.16, 0.84))
  print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))
  
print("\nSyserr of continuum")
for i in range(ncode):
  mean, low, up = np.quantile(sample[:, num_params_var+2*ncode+i], q=(0.5, 0.16, 0.84))
  print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))

print("\nlog10 Error Scale of continuum")
for i in range(ncode):
  mean, low, up = np.quantile(sample[:, num_params_var+3*ncode+i], q=(0.5, 0.16, 0.84))
  print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))

if config["dump"]["fline"] != "":
  print("\nSyserr of line")
  for i in range(ncode):
    mean, low, up = np.quantile(sample[:, num_params_var+4*ncode+i], q=(0.5, 0.16, 0.84))
    print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))

  print("\nlog10 Error Scale of line")
  for i in range(ncode):
    mean, low, up = np.quantile(sample[:, num_params_var+5*ncode+i], q=(0.5, 0.16, 0.84))
    print(code[i], "%5.3f -%5.3f +%5.3f"%(mean, mean-low, up-mean))
 
plt.rc('text', usetex=True)
plt.rc('font', family="serif", size=18)

#===================================================================
# now plot
#===================================================================
data={}
nax = 1
# first continuum 
cont = np.loadtxt(basename(config["dump"]["fcont"]))
cont_code_org = np.empty(cont.shape[0], dtype="U20")
cont_cali = np.loadtxt(basename(config["dump"]["fcont"])+"_cali", usecols=(0, 1, 2))
cont_code = np.loadtxt(basename(config["dump"]["fcont"])+"_cali", usecols=(3), dtype=str)
data["cont"]=[cont, cont_cali]
cont_full = np.loadtxt("cont_recon.txt")

# create original code of the raw data
i1=0
i2=0
for i in range(ncode):
  i2 = i1 + np.count_nonzero(cont_code==code[i])
  cont_code_org[i1:i2]=code[i]
  i1 = i2

# load index for sorting the data
idx_cont = np.loadtxt("cont_sort_index.txt", dtype=int)

# load line data if included
if config["dump"]["fline"] != "":
  nax+=1
  line = np.loadtxt(basename(config["dump"]["fline"]))
  line_code_org = np.empty(line.shape[0], dtype="U20")
  line_cali = np.loadtxt(basename(config["dump"]["fline"])+"_cali", usecols=(0, 1, 2))
  line_code = np.loadtxt(basename(config["dump"]["fline"])+"_cali", usecols=(3), dtype=str)
  data["line"] = [line, line_cali]
  line_full = np.loadtxt("line_recon.txt")
  idx_line = np.loadtxt("line_sort_index.txt", dtype=int)
  
  i1=0
  i2=0
  for i in range(ncode):
    i2 = i1 + np.count_nonzero(line_code==code[i])
    line_code_org[i1:i2]=code[i]
    i1 = i2

# obtain colors of matplotlib
#cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
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
              ecolor='grey', markeredgecolor=None, elinewidth=1, capsize=0.9,  label=r'${0}~({1})$'.format(code[i], len(idx[0])))
              
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
  ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=1.5)
  ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=d[idx_cont[idx[0]], 2], ls='none', color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=1.5)
  
ax.set_ylabel("Intercalibrated Flux")
ax.set_xlim(xlim[0], xlim[1])
[xt.set_visible(False) for xt in ax.get_xticklabels()]
ax.minorticks_on()

# plot parameter prior
ax = fig.add_axes((0.76, 0.38, 0.2, 0.5))
ax.text(0.3, 0.5, r"$\varphi,~~G, ~~\epsilon, ~~b$", fontsize=15)
for i in range(ncode):
  fstr = r"${0}$".format(code[i])
  ax.text(0.1, 0.45-i*0.04, fstr, fontsize=15)
  fstr=r""
  if np.std(sample[:, num_params_var+i]) == 0.0 :
    fstr = fstr + r"1"
  else:
    fstr = fstr + r"0"
  
  if np.std(sample[:, num_params_var+ncode+i]) == 0.0 :
    fstr = fstr + r"~~~~~1"
  else:
    fstr = fstr + r"~~~~~0"
  
  if np.std(sample[:, num_params_var+2*ncode+i]) == 0.0 :
    fstr = fstr + r"~~~~1"
  else:
    fstr = fstr + r"~~~~0"
  
  if np.std(sample[:, num_params_var+3*ncode+i]) == 0.0 :
    fstr = fstr + r"~~~1"
  else:
    fstr = fstr + r"~~~0"
  
  ax.text(0.3, 0.45-i*0.04, fstr, fontsize=15)

ax.text(0.1, 0.45-ncode*0.04, "0: free, 1: fixed", fontsize=15)
ax.set_axis_off()

# plot residuals
ax = fig.add_axes((0.1, 0.08, 0.66, 0.28))
for i in range(ncode):
 idx = np.where((cont_code == code[i]))
 res = dc[idx[0], 1] - np.interp(dc[idx[0], 0], cont_full[:, 0], cont_full[:, 1])
 ax.errorbar(dc[idx[0], 0], res, yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor=cycle[np.mod(i, len(cycle))], markeredgecolor=None,  elinewidth=1, capsize=1.5,  label=r'${0}$'.format(code[i]), zorder=1)
 
 ax.errorbar(dc[idx[0], 0], res, yerr=d[idx_cont[idx[0]], 2], ls='none', color=cycle[np.mod(i, len(cycle))], \
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
  ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(dc[idx[0], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
             elinewidth=1, capsize=1.5, zorder=1)
  ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(d[idx_cont[idx[0]], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
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
ax.hist((dc[:, 1] - np.interp(dc[:, 0], cont_full[:, 0], cont_full[:, 1]))/dc[:, 2], orientation='horizontal', density=True, bins=10)
y = np.linspace(-4, 4, 100)
x = 1.0/np.sqrt(2.0*np.pi)*np.exp(-0.5*y*y)
ax.plot(x, y)
ax.set_ylim(-4, 4)
#[yt.set_visible(False) for yt in ax.get_yticklabels()]
ax.set_ylabel("Stardarized Residuals")
ax.minorticks_on()

fname = basename(config["dump"]["fcont"])
fname = fname.replace("_", " ")
fig.suptitle(r"\bf {0}".format(fname), x=0.5, y=1.0)
pdf.savefig(fig)
plt.close()

#===================================================================
# then plot line if there is
#===================================================================
if config["dump"]["fline"] != "":
 fig = plt.figure(figsize=(15, 12))
 
 ax = fig.add_axes((0.1, 0.68, 0.66, 0.28))
 key="line"
 d = data[key][0]
 dc = data[key][1]
 for i in range(ncode):
  idx = np.where((line_code_org == code[i]))
  ax.errorbar(d[idx[0], 0], d[idx[0], 1], yerr=d[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor='grey', markeredgecolor=None, elinewidth=1, capsize=0.9,  label=r'${0}~({1})$'.format(code[i], len(idx[0])))
 
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
  ax.errorbar(dc[idx[0], 0], dc[idx[0], 1], yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor='grey', markeredgecolor=None,  elinewidth=1, capsize=0.9, label=r'${0}$'.format(code[i]))
 
 ax.set_ylabel("Intercalibrated Flux")
 ax.set_xlim(xlim[0], xlim[1])
 ax.minorticks_on()
 [xt.set_visible(False) for xt in ax.get_xticklabels()]

 # plot parameter prior
 ax = fig.add_axes((0.76, 0.38, 0.2, 0.5))
 ax.text(0.3, 0.5, r"$\varphi,~~G, ~~\epsilon, ~~b$", fontsize=15)
 for i in range(ncode):
   fstr = r"${0}$".format(code[i])
   ax.text(0.1, 0.45-i*0.04, fstr, fontsize=15)
   fstr=r""
   if np.std(sample[:, num_params_var+i]) == 0.0 :
     fstr = fstr + r"1"
   else:
     fstr = fstr + r"0"
   
   if np.std(sample[:, num_params_var+ncode+i]) == 0.0 :
     fstr = fstr + r"~~~~~1"
   else:
     fstr = fstr + r"~~~~~0"
   
   if np.std(sample[:, num_params_var+2*ncode+i]) == 0.0 :
     fstr = fstr + r"~~~~1"
   else:
     fstr = fstr + r"~~~~0"
   
   if np.std(sample[:, num_params_var+3*ncode+i]) == 0.0 :
     fstr = fstr + r"~~~1"
   else:
     fstr = fstr + r"~~~0"
   
   ax.text(0.3, 0.45-i*0.04, fstr, fontsize=15)

 ax.text(0.1, 0.45-ncode*0.04, "0: free, 1: fixed", fontsize=15)
 ax.set_axis_off()


 ax = fig.add_axes((0.1, 0.08, 0.66, 0.28))
 for i in range(ncode):
   idx = np.where((line_code == code[i]))
   res = dc[idx[0], 1] - np.interp(dc[idx[0], 0], line_full[:, 0], line_full[:, 1])
   ax.errorbar(dc[idx[0], 0], res, yerr=dc[idx[0], 2], ls='none', marker='o', markersize=3, color=cycle[np.mod(i, len(cycle), dtype=int)], \
              ecolor=cycle[np.mod(i, len(cycle), dtype=int)], markeredgecolor=None,  elinewidth=1, capsize=0.9, label=r'${0}$'.format(code[i]), zorder=0)
   
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
   ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(dc[idx[0], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
              elinewidth=1, capsize=1.5, zorder=1)
   ax.errorbar(xlim[1]-(xlim[1]-xlim[0])/(ncode+4) * (i+2), 0.0, yerr=np.mean(d[idx_line[idx[0]], 2]), color=cycle[np.mod(i, len(cycle), dtype=int)],\
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
 ax.hist((dc[:, 1] - np.interp(dc[:, 0], line_full[:, 0], line_full[:, 1]))/dc[:, 2], orientation='horizontal', density=True, bins=10)
 y = np.linspace(-4, 4, 100)
 x = 1.0/np.sqrt(2.0*np.pi)*np.exp(-0.5*y*y)
 ax.plot(x, y)
 ax.set_ylim(-4, 4)
 
 #[yt.set_visible(False) for yt in ax.get_yticklabels()]
 ax.set_ylabel("Stardarized Residuals")
 ax.minorticks_on()
 
 fname = basename(config["dump"]["fline"])
 fname = fname.replace("_", " ")
 fig.suptitle(r"\bf {0}".format(fname), x=0.5, y=1.0)

 pdf.savefig(fig)
 plt.close()

#pdf.close()
#exit()

#===================================================================
# now plot histograms
#===================================================================
if int(config["dump"]["fixed_scale"]) == 1:
  fig = corner.corner(sample[:, num_params_var+ncode+1:num_params_var+2*ncode], smooth=True, smooth1d = True, \
        levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
  
  ax = fig.get_axes()
  for i in range(ncode-1):
    xlim = ax[i*(ncode-1)+i].get_xlim()
    ylim = ax[i*(ncode-1)+i].get_ylim()
    ax[i*(ncode-1)+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i+1]))
  fig.suptitle(r"\bf Shift", fontsize=20)
  pdf.savefig(fig)
  plt.close()
 
elif int(config["dump"]["fixed_shift"]) == 1:
  fig = corner.corner(sample[:, num_params_var+1:num_params_var+ncode], smooth=True, smooth1d = True,  \
        levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
  ax = fig.get_axes()
  for i in range(ncode-1):
    xlim = ax[i*(ncode-1)+i].get_xlim()
    ylim = ax[i*(ncode-1)+i].get_ylim()
    ax[i*(ncode-1)+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i+1]))
  fig.suptitle(r"\bf Scale", fontsize=20)
  pdf.savefig(fig)
  plt.close()
else:
  for i in range(1, ncode):
    range_min = np.min(sample[:, [num_params_var+i,num_params_var+i+ncode]], axis=0)
    range_max = np.max(sample[:, [num_params_var+i,num_params_var+i+ncode]], axis=0)
    span = range_max - range_min
    range_interval = [[range_min[i]-0.3*span[i], range_max[i]+0.3*span[i]] for i in range(2)]
    fig = corner.corner(sample[:, [num_params_var+i,num_params_var+i+ncode]], smooth=True, smooth1d = True, labels=[r"$\log\varphi$", r"$G$"], 
          range=range_interval, \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
    
    ax = fig.get_axes()
    ax[1].text(0.0, 0.5, r"\bf Scale \& Shift", fontsize=15)
    ax[1].text(0.0, 0.65, r"\bf"+r'$\bf {0}$'.format(code[i]), fontsize=15)
    #fig.suptitle(r"\bf Scale \& Shift  "+r'${0}$'.format(code[i][3:-4]), fontsize=20)
    pdf.savefig(fig)
    plt.close()

if int(config["dump"]["fixed_syserr"]) == 0:
  fig = corner.corner(sample[:, num_params_var+2*ncode:num_params_var+3*ncode], smooth=True, smooth1d = True, \
       levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
  ax = fig.get_axes()
  for i in range(ncode):
    xlim = ax[i*ncode+i].get_xlim()
    ylim = ax[i*ncode+i].get_ylim()
    ax[i*ncode+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i]))
  fig.suptitle(r"\bf Systematic Error (Continuum)", fontsize=20)
  pdf.savefig(fig)
  plt.close()


if int(config["dump"]["fixed_error_scale"]) == 0:
  fig = corner.corner(sample[:, num_params_var+3*ncode:num_params_var+4*ncode], smooth=True, smooth1d = True, show_titles=True, title_fmt=".3f")
  ax = fig.get_axes()
  for i in range(ncode):
    xlim = ax[i*ncode+i].get_xlim()
    ylim = ax[i*ncode+i].get_ylim()
    ax[i*ncode+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i]))
  fig.suptitle(r"\bf Error Scale", fontsize=20)
  pdf.savefig(fig)
  plt.close()

if int(config["dump"]["fixed_syserr"]) == 0 and int(config["dump"]["fixed_error_scale"]) == 0:

  for i in range(ncode):
    fig = corner.corner(sample[:, [num_params_var+2*ncode+i,num_params_var+2*ncode+i+ncode]], smooth=True, smooth1d = True, labels=[r"$\epsilon$", r"$b$"], 
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
    
    ax = fig.get_axes()
    ax[1].text(0.0, 0.5, r"\bf Syserr \& Error Scale", fontsize=15)
    ax[1].text(0.0, 0.6, r"\bf"+r'${0}$'.format(code[i]), fontsize=15)
    #fig.suptitle(r"\bf Syserr \& Error Scale  "+code[i][3:-4], fontsize=20)
    pdf.savefig(fig)
    plt.close()

if config["dump"]["fline"] != "":
  if int(config["dump"]["fixed_syserr"]) == 0:
    fig = corner.corner(sample[:, num_params_var+4*ncode:num_params_var+5*ncode], smooth=True, smooth1d = True, \
          levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
    ax = fig.get_axes()
    for i in range(ncode):
      xlim = ax[i*ncode+i].get_xlim()
      ylim = ax[i*ncode+i].get_ylim()
      ax[i*ncode+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i]))
    fig.suptitle(r"\bf Systematic Error (Line)", fontsize=20)
    pdf.savefig(fig)
    plt.close()

  if int(config["dump"]["fixed_error_scale"]) == 0:
    fig = corner.corner(sample[:, num_params_var+5*ncode:num_params_var+6*ncode], smooth=True, smooth1d = True,\
      levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
    ax = fig.get_axes()
    for i in range(ncode):
      xlim = ax[i*ncode+i].get_xlim()
      ylim = ax[i*ncode+i].get_ylim()
      ax[i*ncode+i].text(xlim[1]-0.2*(xlim[1]-xlim[0]), ylim[1] - 0.2*(ylim[1]-ylim[0]), r'$\bf {0}$'.format(code[i]))
    fig.suptitle(r"\bf Error Scale (Line)", fontsize=20)
    pdf.savefig(fig)
    plt.close()
    
  if int(config["dump"]["fixed_syserr"]) == 0 and int(config["dump"]["fixed_error_scale"]) == 0:

    for i in range(ncode):
      fig = corner.corner(sample[:, [num_params_var+4*ncode+i,num_params_var+4*ncode+i+ncode]], smooth=True, smooth1d = True, \
          labels=[r"$\epsilon$", r"$b$"], levels=1.0-np.exp(-0.5*np.arange(1.0, 3.1, 1.0)**2), show_titles=True, title_fmt=".3f")
    
      ax = fig.get_axes()
      ax[1].text(0.0, 0.5, r"\bf Syserr \& Error Scale", fontsize=15)
      ax[1].text(0.0, 0.65, r"\bf"+r'$\bf {0}$'.format(code[i]), fontsize=15)
      pdf.savefig(fig)
      plt.close()

  
pdf.close()

