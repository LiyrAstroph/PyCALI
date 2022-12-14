#
# generate formatted input files for PyCALI.
#
import pathlib
import numpy as np
import pandas as pd

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

def format(fname, data):
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
      
    fp.write("# %s %d\n"%(key, len(d[:, 0])))
    for i in range(len(d[:, 0])):
      fp.write("%16.10e %e %e\n"%(d[i, 0], d[i, 1], d[i, 2]))
    
  fp.close()

def convert_asassn(datafile, unit=3.92e-9, time_start=0.0, rebin=False, errlimit=0.1, diffcamera=False):
  """
  convert asassn  datafile into flux

  unit=3.92e-9 is the V-band zero flux point

  """  
  if not isinstance(datafile, str):
    raise ValueError("Input datafile is not a string!")
  
  path = pathlib.Path(datafile)
  if path.suffix != ".csv":
    raise ValueError("fname is not a csv file.")
  
  asas_all = np.loadtxt(datafile, delimiter=',', usecols=(0, 5, 6), skiprows=1)
  band = np.loadtxt(datafile, delimiter=',', usecols=(2, 9), skiprows=1, dtype=str)
  # first sort data
  arg = np.argsort(asas_all[:, 0])
  asas_all = asas_all[arg, :]
  band = band[arg]
  
  # check error limit 
  idx = np.where(asas_all[:, 2]<=errlimit)
  asas_all = asas_all[idx[0], :]
  band = band[idx[0]]
  
  # convert to flux
  asas_all[:, 0] -= time_start
  asas_all[:, 1] = 10.0**(-asas_all[:, 1]/2.5) / unit
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
          
        if rebin == True:
          tc, yc, yerrc = data_rebin(asas_all[idx[0], 0], asas_all[idx[0], 1], asas_all[idx[0], 2], 1)
          asas["asas_"+c+f] = np.stack((tc, yc, yerrc), axis=-1)
        else:
          asas["asas_"+c+f] = asas_all[idx[0], :]
    else:
      idx = np.where(band[:, 1]==f)
      if len(idx[0]) == 0:
        continue
        
      if rebin == True:
        tc, yc, yerrc = data_rebin(asas_all[idx[0], 0], asas_all[idx[0], 1], asas_all[idx[0], 2], 1)
        asas["asas_"+f] = np.stack((tc, yc, yerrc), axis=-1)
      else:
        asas["asas_"+f] = asas_all[idx[0], :]
  
  return asas


def convert_ztf(datafile, unit=3.92e-9, time_start=0.0, rebin=False, errlimit=0.1):
  """
  convert ZTF datafile into flux

  unit=3.92e-9 is the V-band zero flux point

  """  
  if not isinstance(datafile, str):
    raise ValueError("Input datafile is not a string!")

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
  mag = 10.0**(-mag/2.5) / unit
  err = mag * err/2.5 * np.log(10.0)

  ztf = {}
  for f in filt_unique:
    idx = np.where(filt==f)
    if len(idx[0]) == 0:
      continue
        
    key = "ztf_"+f
    if rebin == True:
      tc, yc, yerrc = data_rebin(jd[idx[0]], mag[idx[0]], err[idx[0]], 1)
      ztf[key] = np.stack((tc, yc, yerrc), axis=-1)
    else:
      ztf[key] = np.stack((jd[idx[0]], mag[idx[0]],err[idx[0]]), axis=-1)

  return ztf

if __name__ == "__main__":
  ztf = convert_ztf("ZTF.csv", rebin=True, errlimit=0.08)
  asassn = convert_asassn("asas.csv", rebin=True, errlimit=0.08)

  data = ztf | asassn

  format("test.txt", data)
