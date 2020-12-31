#
# generate mock data for pyCali.
#
import numpy as np
from numpy import fft 
import pycali
import matplotlib.pyplot as plt 
import copy 

def convolve_fft(con, resp):
  """
  convolution continuum with a response function using FFT
  """
  resp_pad = np.zeros(con.shape[0])
  resp_pad[:resp.shape[0]] = resp
  con_fft = fft.rfft(con)
  resp_fft = fft.rfft(resp_pad)
  conv_fft = con_fft * resp_fft 
  conv = fft.irfft(conv_fft, n = con.shape[0])
  return conv

def generate_mock_data():
  """
  generate mock data
  """
  # DRW parameters
  sigma = 0.3
  tau = 50.0
  
  # time nodes for continuum
  tg = np.linspace(-200.0, 600.0, 2000)
  tmx, tmy = np.meshgrid(tg, tg)
  # covariance
  Cov = sigma * sigma * np.exp(-np.abs(tmx-tmy)/tau)
  Mat = np.linalg.cholesky(Cov)
  u = np.random.randn(tg.shape[0])
  fs = np.matmul(Mat, u) + 1.0
  # errors, around 0.015
  fe = np.random.randn(tg.shape[0])*0.005+0.015
  con = np.stack((tg, fs, fe), axis=-1)
  
  # now emission line 
  # first transfer function, Gaussian
  dt = con[1, 0] - con[0, 0]
  ntau = 200
  resp = np.zeros(ntau)
  tau = np.array(np.arange(ntau))*dt
  resp[:] = np.exp(-0.5 * (tau-30.0)**2/(10.0)**2)
  resp[:] /= np.sum(resp[:]) * dt
  
  # get emission line
  nline = 1000
  conv = convolve_fft(con[:, 1], resp) * dt
  line = np.zeros((nline, 3))
  line[:, 0] = np.linspace(0.0, 360.0, nline)
  line[:, 2] = np.random.randn(line.shape[0])*0.005+0.015
  line[:, 1] = np.interp(line[:, 0], con[:, 0], conv) 
  

  # section between 0-360day
  idx = np.where((con[:, 0]>=0.0) & (con[:, 0]<=360.0))
  con = con[idx[0], :]
  
  # code, N, scale, shift etc
  codes=["A", "B", "C", "D", "E"]
  num_cont = [150, 120, 100, 82, 180]
  num_line = [160, 90, 120, 80, 100]
  scale = [1.0, 0.9, 1.1, 0.94, 1.05]
  shift = [0.0, -0.2, 0.15, 0.05, -0.11]
  syserr_cont = [0.0, 0.01, 0.09, 0.015, 0.02]
  syserr_line = [0.0, 0.01, 0.09, 0.015, 0.02]
  error_scale_cont = np.array([1.0, 0.7, 0.9, 1.2, 0.8])*0.0 + 1.0
  error_scale_line = np.array([1.0, 1.3, 0.8, 1.3, 1.1])*0.0 + 1.0
  
  print("code:", codes)
  print("scale:", scale)
  print("shift:", shift)
  print("syserr cont:", syserr_cont)
  print("syserr line:", syserr_line)
  print("error scale cont:", error_scale_cont)
  print("error scale line:", error_scale_line)
  
  fig = plt.figure()
  ax1 = fig.add_subplot(211)
  ax2 = fig.add_subplot(212)

  # save full continuum 
  np.savetxt("data/sim_cont_full.txt", con, fmt="%15.5f")
  fp=open("data/sim_cont.txt", "w")
  ax1.errorbar(con[:, 0], con[:, 1], yerr=con[:, 2], ls='none')
  for i in range(len(codes)):
    idx = np.unique(np.random.randint(con.shape[0],size=num_cont[i]))
    con_set = con[idx, :]
    # flux + noise + systematic error
    con_set[:, 1] = (con_set[:, 1] + shift[i])/scale[i] + np.random.randn(con_set.shape[0]) * con_set[:, 2] \
                  + np.random.randn(con_set.shape[0])*syserr_cont[i]
    con_set[:, 2] = con_set[:, 2]/error_scale_cont[i]
    ax1.errorbar(con_set[:, 0], con_set[:, 1], yerr=con_set[:, 2], ls='none')
    
    fp.write("# %s %d\n"%(codes[i], con_set.shape[0]))
    np.savetxt(fp, con_set, fmt="%15.5f")

  fp.close()
  ax1.set_ylabel('Continuum')
  ax1.minorticks_on()

  np.savetxt("data/sim_line_full.txt", line, fmt="%15.5f")
  fp=open("data/sim_line.txt", "w")
  ax2.errorbar(line[:, 0], line[:, 1], yerr=line[:, 2], ls='none')
  for i in range(len(codes)):
    idx = np.unique(np.random.randint(line.shape[0],size=num_line[i]))
    line_set = line[idx, :]
    line_set[:, 1] = (line_set[:, 1])/scale[i] + np.random.randn(line_set.shape[0]) * line_set[:, 2] \
                   + np.random.randn(line_set.shape[0])*syserr_line[i]
    line_set[:, 2] = line_set[:, 2]/error_scale_line[i]
    ax2.errorbar(line_set[:, 0], line_set[:, 1], yerr=line_set[:, 2], ls='none')
    
    fp.write("# %s %d\n"%(codes[i], line_set.shape[0]))
    np.savetxt(fp, line_set, fmt="%15.5f")
  
  fp.close()
  ax2.set_xlabel('Time (day)')
  ax2.set_ylabel('Line')
  ax2.minorticks_on()
  fig.suptitle("Mock Data")
  plt.show()


if __name__ == "__main__":
  generate_mock_data()