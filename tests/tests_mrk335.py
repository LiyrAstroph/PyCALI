#!/usr/bin/env python
# coding: utf-8

# <font size=4>**Here is a detailed tutorial to show how to use PyCALI in Python by taking the data of Mrk 335 from ASAS-SN and ZTF as an example.**</font>

# <font size=4>**First download the data of Mrk 335 in a CSV format from ASAS-SN website (https://asas-sn.osu.edu/) and ZTF website (ZTF g band; https://irsa.ipac.caltech.edu/Missions/ztf.html). Suppose the filenames are Mrk335_asas.csv and Mrk335_ztf.csv, respectively.** </font>
# 
# <font size=4>**Load these data and generate a formatted input file for PyCALI.** </font>

# In[14]:


import numpy as np
import pycali
import matplotlib.pyplot as plt


# In[3]:


ztf = pycali.convert_ztf("Mrk335_ztf.csv", rebin=True, errlimit=0.079, unit=3.92e-9)
# rebin:  whether rebin the points within one day
# errlimit: discard these points with errors larger than this limit
# unit is the zero-magnitude flux density
# return a dict
    
asas = pycali.convert_asassn("Mrk335_asas.csv", rebin=True, errlimit=0.079, diffcamera=False, unit=3.92e-9)
# diffcamera: whether treat different cameras as different datasets
    
data = ztf | asas  # combine the two dicts
    
# write to a file named "Mrk335.txt"
pycali.format("Mrk335.txt", data)

# if only use data in a time range, use
# pycali.format("Mrk335.txt", data, trange=(2458200, 2470000))


# In[11]:


print("asas:", asas.keys())
print("ztf:", ztf.keys())
print("data:", data.keys())


# <font size=4>**Now the input file has been created. Next run PyCALI to do intercalibration.**</font>
# 
# <font size=4>**This will take a while to run ..., better to directly run the python script "tests_mrk335.py", which shows the running status in the terminal.**</font>

# In[12]:


# setup configurations
cfg = pycali.Config()
cfg.setup(
      fcont="./Mrk335.txt",     # fcont is a string 
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
pycali.plot_results(cfg)
    
# a simple plot 
pycali.simple_plot(cfg)


# <font size=4>**Take a look at the intercalibrated data**</font>

# In[18]:


data_cali = np.loadtxt("Mrk335.txt_cali", usecols=(0, 1, 2))
code = np.loadtxt("Mrk335.txt_cali", usecols=(3), dtype=str)
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)
for c in np.unique(code):
    idx = np.where(code == c)[0]
    ax.errorbar(data_cali[idx, 0],  data_cali[idx, 1], yerr=data_cali[idx, 2], ls='none', marker='o', markersize=3, label=c)

ax.legend()
ax.set_title("Intercalibrated data")
plt.show()


# <font size=4>**There appears a number of outliers. We can discard these outliers as follows.**</font>

# In[21]:


# load data
data = pycali.load_pycali_data("Mrk335.txt")

# load intercalibrated data and ancillary files
cali = np.loadtxt("Mrk335.txt_cali", usecols=(0, 1, 2))
code = np.loadtxt("Mrk335.txt_cali", usecols=(3), dtype=str)
recon = np.loadtxt("Mrk335.txt_recon")
intp = np.interp(cali[:, 0], recon[:, 0], recon[:, 1])
err = np.interp(cali[:, 0], recon[:, 0], recon[:, 2])

# residuals between the calibrated data and reconstruction with a DRW process
res = (cali[:, 1]-intp)/err

# now delete bad points with residual > 5 sigma
data_new = {}
code_uni = np.unique(code)
for c in code_uni:
  idx = np.where((code == c))[0]
  res_code = res[idx]
  idx = np.where(np.abs(res_code)>5)[0]
  data_new[c] = np.delete(data[c], idx, 0)

pycali.format("Mrk335_new.txt", data_new)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)
plt.plot(cali[:, 0], res, ls='none', marker='o')
plt.axhline(y=5, ls='--')
plt.axhline(y=-5, ls='--')
ax.set_ylabel("Res")
ax.set_title("Standarized residuals")
plt.show()


# <font size=4>**Now redo the intercalibration on new data.**</font>

# In[22]:


# setup configurations
cfg = pycali.Config()
cfg.setup(
      fcont="./Mrk335_new.txt",     # fcont is a string 
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
pycali.plot_results(cfg)
    
# a simple plot 
pycali.simple_plot(cfg)


# <font size=4>**Take a look at the newly intercalibrated data**</font>

# In[24]:


data_cali_new = np.loadtxt("Mrk335_new.txt_cali", usecols=(0, 1, 2))
code = np.loadtxt("Mrk335_new.txt_cali", usecols=(3), dtype=str)
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)
for c in np.unique(code):
    idx = np.where(code == c)[0]
    ax.errorbar(data_cali_new[idx, 0],  data_cali_new[idx, 1], yerr=data_cali_new[idx, 2], ls='none', marker='o', markersize=3, label=c)

ax.legend()
ax.set_title("Intercalibrated data")
plt.show()


# In[ ]:




