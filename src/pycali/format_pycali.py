#
# generate formatted input files for PyCALI.
#
import pathlib
import numpy as np
import pandas as pd

__all__=["data_rebin", "format", "convert_ztf", "convert_asassn", "convert_mydata", 
         "remove_outliers", "load_pycali_data", "load_pycali_data_flag"]

def data_rebin(x, y, ye, tb):
  i = 0
  ic = 0
  xc = np.zeros(len(x), dtype=np.double)
  yc = np.zeros(len(x), dtype=np.double)
  yerr = np.zeros(len(x), dtype=np.double)
  
  while i < len(x):
    tmean = 0.0
    fmean = 0.0
    norm = 0.0
    nc = 0
    for j in range(i, len(x)):
      if x[j] < x[i] + tb:
        ye2 = ye[j]*ye[j]
        tmean += x[j]/(ye2)
        fmean += y[j]/(ye2)
        norm += 1.0/ye2
        nc += 1
      
    jp = nc + i - 1
           
    tmean /= norm
    fmean /= norm
    error = np.sqrt(1.0/norm)
      
    sig = error
    
    xc[ic] = tmean
    yc[ic] = fmean
    yerr[ic] = sig
    
    i = jp  + 1
    ic += 1
  
  return xc[:ic], yc[:ic], yerr[:ic]

def format(fname, data, trange=None, unit=1.0, time_start=0.0):
  """
  generate PyCALI formatted file for data
  """
  if not isinstance(fname, str):
    raise ValueError("Input fname is not a string!")
  
  if not isinstance(data, dict):
    raise ValueError("Input data is not dict!")
  
  fp = open(fname, "w")
  
  for key in data.keys():
    d = data[key]   
    
    if len(d[:, 0]) == 0:
      continue

    if trange == None:  
      fp.write("# %s %d\n"%(key, len(d[:, 0])))
      for i in range(len(d[:, 0])):
        fp.write("%16.10e %e %e\n"%(d[i, 0]-time_start, d[i, 1]/unit, d[i, 2]/unit))
    else:
      idx = np.where((d[:, 0]>=trange[0])&(d[:, 0]<=trange[1]))[0]
      if len(idx) == 0:
        continue
      
      fp.write("# %s %d\n"%(key, len(idx)))
      for i in idx:
        fp.write("%16.10e %e %e\n"%(d[i, 0]-time_start, d[i, 1]/unit, d[i, 2]/unit))
    
  fp.close()

def convert_asassn(datafile, useflux=False, zeropoint=3.92e-9, time_start=0.0, rebin=None, errlimit=0.1, diffcamera=False, keylabel=""):
  """
  convert asassn  datafile into flux

  unit=3.92e-9 is the V-band zero flux point

  """  
  if not isinstance(datafile, str):
    raise ValueError("Input datafile is not a string!")
  
  path = pathlib.Path(datafile)
  if path.suffix != ".csv":
    raise ValueError("fname is not a csv file.")
  
  # binning interval
  is_rebin = False
  rebin_interval = 1
  if rebin is None:
    is_rebin = False 
  elif type(rebin) == "bool": # using 1 day
    is_rebin = rebin 
    rebin_interval = 1
  else:  # otherwise, using input value
    is_rebin = True 
    rebin_interval = float(rebin)

  data = pd.read_csv(datafile, comment="#")
  
  band = np.column_stack((np.array(data["Camera"], dtype=str), np.array(data["Filter"], dtype=str)))
  #band = np.genfromtxt(datafile, delimiter=',', usecols=(2, 9), skip_header=1, dtype=str)

  if useflux == False:
    #asas_all = np.genfromtxt(datafile, delimiter=',', usecols=(0, 5, 6), skip_header=1)
    if "Flux" in data.keys():  # Sky Patrol V2
      asas_all = np.column_stack((data["JD"], data["Mag"], data["Mag Error"]))
    else:
      # remove bad points 
      idx = np.where(data["mag_err"]!=99.99)
      asas_all = np.column_stack((data["HJD"][idx[0]], np.array(data["mag"][idx[0]], dtype=float), 
                                  data["mag_err"][idx[0]]))
      band = band[idx[0], :]
      
  else:
    #asas_all = np.genfromtxt(datafile, delimiter=',', usecols=(0, 7, 8), skip_header=1)
    if "Mag" in data.keys(): # Sky Patrol V2
      asas_all = np.column_stack((data["JD"], data["Flux"], data["Flux Error"]))
    else:
      asas_all = np.column_stack((data["HJD"], data["flux(mJy)"], data["flux_err"]))
      
    # mJy to erg/s/cm^2/A
    # V band
    idx = np.where(band[:, 1] == "V")
    asas_all[idx[0], 1:] *= 1.0e-26 * 3e10/(5500*1.0e-8)**2 /1.0e8
    # g band
    idx = np.where(band[:, 1] == "g")
    asas_all[idx[0], 1:] *= 1.0e-26 * 3e10/(5200*1.0e-8)**2 /1.0e8

  # remove bad values
  idx = np.logical_not(np.isnan(asas_all[:, 1]))
  asas_all = asas_all[idx, :]
  band = band[idx]

  # first sort data
  arg = np.argsort(asas_all[:, 0])
  asas_all = asas_all[arg, :]
  band = band[arg]
  
  # check error limit 
  if useflux == False: # magnitude
    idx = np.where(asas_all[:, 2]<=errlimit)
  else:  # flux
    idx =  np.where(asas_all[:, 2]/asas_all[:, 1]<=errlimit)

  asas_all = asas_all[idx[0], :]
  band = band[idx[0]]
  
  # convert to flux
  if useflux == False:
    asas_all[:, 0] -= time_start
    asas_all[:, 1] = 10.0**(-asas_all[:, 1]/2.5) * zeropoint
    asas_all[:, 2] = asas_all[:, 1] * asas_all[:, 2]/2.5 * np.log(10.0)
  
  camera = np.unique(band[:, 0])
  filt = np.unique(band[:, 1])

  asas = {}
  for f in filt:
    if diffcamera == True:
      for c in camera:
        idx = np.where((band[:, 0]==c)&(band[:, 1]==f))
        if len(idx[0]) == 0:
          continue
          
        if is_rebin == True:
          tc, yc, yerrc = data_rebin(asas_all[idx[0], 0], asas_all[idx[0], 1], asas_all[idx[0], 2], rebin_interval)
          asas[keylabel+"asas_"+c+f] = np.stack((tc, yc, yerrc), axis=-1)
        else:
          asas[keylabel+"asas_"+c+f] = asas_all[idx[0], :]
    else:
      idx = np.where(band[:, 1]==f)
      if len(idx[0]) == 0:
        continue
        
      if is_rebin == True:
        tc, yc, yerrc = data_rebin(asas_all[idx[0], 0], asas_all[idx[0], 1], asas_all[idx[0], 2], rebin_interval)
        asas[keylabel+"asas_"+f] = np.stack((tc, yc, yerrc), axis=-1)
      else:
        asas[keylabel+"asas_"+f] = asas_all[idx[0], :]
  
  return asas


def convert_ztf(datafile, zeropoint=3.92e-9, time_start=0.0, rebin=None, errlimit=0.1, keylabel=""):
  """
  convert ZTF datafile into flux

  unit=3.92e-9 is the V-band zero flux point

  """  
  if not isinstance(datafile, str):
    raise ValueError("Input datafile is not a string!")
  
  # binning interval
  is_rebin = False
  rebin_interval = 1
  if rebin is None:
    is_rebin = False 
  elif type(rebin) == "bool": # using 1 day
    is_rebin = rebin 
    rebin_interval = 1
  else:  # otherwise, using input value
    is_rebin = True 
    rebin_interval = float(rebin)

  path = pathlib.Path(datafile)
  if path.suffix != ".csv":
    raise ValueError("fname is not a csv file.")
  
  data = pd.read_csv(datafile)
  
  flag = data['catflags'].values
  
  # check ztf flag
  filt = data['filtercode'][flag==0].values
  jd = data['hjd'][flag==0].values
  mag = data['mag'][flag==0].values
  err = data['magerr'][flag==0].values
  
  # check error limit
  idx = np.where(err<=errlimit)
  filt = filt[idx[0]]
  jd = jd[idx[0]]
  mag = mag[idx[0]]
  err = err[idx[0]]

  # first sort data 
  arg = np.argsort(jd)
  filt = filt[arg]
  jd = jd[arg]
  mag = mag[arg]
  err = err[arg]

  filt_unique = np.unique(filt)

  # convert to flux
  jd -= time_start 
  mag = 10.0**(-mag/2.5) * zeropoint
  err = mag * err/2.5 * np.log(10.0)

  ztf = {}
  for f in filt_unique:
    idx = np.where(filt==f)
    if len(idx[0]) == 0:
      continue
        
    key = keylabel+"ztf_"+f
    if is_rebin == True:
      tc, yc, yerrc = data_rebin(jd[idx[0]], mag[idx[0]], err[idx[0]], rebin_interval)
      ztf[key] = np.stack((tc, yc, yerrc), axis=-1)
    else:
      ztf[key] = np.stack((jd[idx[0]], mag[idx[0]],err[idx[0]]), axis=-1)

  return ztf

def convert_mydata(fname, keylabel=""):
  """
  convert to dict from fname
  """

  if keylabel == "":
    keylabel = "mydata"
  
  data = np.loadtxt(fname, usecols=(0, 1, 2))

  data_dict = {keylabel:data}

  return data_dict

def remove_outliers(fname, dev=5, doplot=False):
  """
  remove outliers with a deviation of (dev) sigma from the reconstruction

  presume that the file "fname" is in PyCALI format, the previously intercalibrated file is "fname_cali"
               and the reconstruction file is "fname_recon".
  """
  
  if fname == None:
    raise ValueError("need to input a file name!")

  # load data
  data = load_pycali_data(fname)

  # load intercalibrated data and ancillary files
  try:
    cali = np.loadtxt(fname+"_cali", usecols=(0, 1, 2))
    code = np.loadtxt(fname+"_cali", usecols=(3), dtype=str)
  except:
    raise IOError(fname+"_cali error!")
  
  try:
    recon = np.loadtxt(fname+"_recon")
  except:
    raise IOError(fname+"_recon error!")
  
  intp = np.interp(cali[:, 0], recon[:, 0], recon[:, 1])
  err = np.interp(cali[:, 0], recon[:, 0], recon[:, 2])

  # residuals between the calibrated data and reconstruction with a DRW process
  res = (cali[:, 1]-intp)/err

  # now delete bad points with residual > dev sigma
  data_new = {}
  # note here use data.keys() to retain the order of codes
  for c in data.keys():
      idx = np.where((code == c))[0]
      res_code = res[idx]
      idx = np.where(np.abs(res_code)>dev)[0]
      data_new[c] = np.delete(data[c], idx, 0)
  
  path = pathlib.Path(fname)
  fname_new = str(path.parent.joinpath(path.stem+"_new.txt"))
  format(fname_new, data_new)

  # do plotting
  if doplot:
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)
    plt.plot(cali[:, 0], res, ls='none', marker='o')
    plt.axhline(y=dev, ls='--', color='k')
    plt.axhline(y=-dev, ls='--', color='k')
    ax.set_ylabel("Res")
    ax.set_title("Standarized residuals")
    plt.show()


def load_pycali_data(fname):
  """
  load data from a formatted file
  return a dict
  """
  data = {}
  block = []
  nc = 0
  code = ""
  fp = open(fname)
  for line in fp:
    if line[0] == "#":
      code = line[1:].split()[0]
      nc = int(line[1:].split()[1])
      block.clear()
    else:
      block.append(np.array(line.split(), dtype=float))

    if len(block) == nc:
      data[code] = np.array(block)

  fp.close()
  return data

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

if __name__ == "__main__":
  ztf = convert_ztf("ZTF.csv", rebin=True, errlimit=0.08)
  asassn = convert_asassn("asas.csv", rebin=True, errlimit=0.08)

  data = ztf | asassn

  format("test.txt", data)
