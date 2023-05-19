**************
Detailed Usage
**************

C/C++ executable binary: cali
-----------------------------

.. code-block:: bash
  
  ./cali param.txt 

Here **param.txt** is the parameter file passed to **cali**. 
The parameter file looks like::

  # parameter file
  # lines beginning with '#' are regarded as comments and are neglected

  FileCont      data/ngc5548_cont.txt
  #FileLine      data/ngc5548_line.txt
  
  #NMcmc         10000
  
  #PTol         0.1
  
  #FixedScale     0    
  #FixedShift     0
  
  #ScaleRangeLow  0.5    
  #ScaleRangeUp   1.5
  
  #ShiftRangeLow  -1.0
  #ShiftRangeUp   1.0
  
  #FixedSyserr     1
  #FixedErrorScale  1
  
  #SyserrRangeLow  0.0  
  #SyserrRangeUp   0.1
  
  #ErrscaleRangeLow  0.1
  #ErrscaleRangeUp   2.0

  #FlagNorm  0
  
  #SigmaRangeLow  1.0e-4
  #SigmaRangeUp  1.0

  #TauRangeLow  1.0
  #TauRangeUp  1.0e4

  #FixedCodes  1,3


In the parameter file, except for the option **FileCont**, all the rest options are optional. If they are not specified, 
**cali** will use the default values as shown above. The meaning of the above options are 

+------------------+-----------------------+---------+--------------------------------------+
| Option           |        Example        | Type    |              Meaning                 |
+==================+=======================+=========+======================================+
| FileCont         | data/ngc5548_cont.txt |         |file name of continuum light curve    |
+------------------+-----------------------+---------+--------------------------------------+
| FileLine         | data/ngc5548_line.txt |optional |file name of line light curves,       |
|                  |                       |         |                                      |
|                  |                       |         |can be multiple lines, use comma (,)  |
|                  |                       |         |to separate, e.g.,                    |
|                  |                       |         |                                      |
|                  |                       |         |**data/ngc5548_line1.txt,**           |
|                  |                       |         |**data/ngc5548_line2.txt**            |
+------------------+-----------------------+---------+--------------------------------------+
| NMcmc            | 10000                 |optional |number of mcmc steps                  |
+------------------+-----------------------+---------+--------------------------------------+
| PTol             | 0.1                   |optional |tolerance of log likelihood in        |
|                  |                       |         |MCMC samling                          |
+------------------+-----------------------+---------+--------------------------------------+
| FixedScale       | 0                     |optional |1: fix scale factor; 0: not           |
+------------------+-----------------------+---------+--------------------------------------+
| FixedShift       | 0                     |optional |1: fix shift factor; 0: not           |
+------------------+-----------------------+---------+--------------------------------------+
| ScaleRangeLow    | 0.5                   |optional |lower limit of scale factor           |
+------------------+-----------------------+---------+--------------------------------------+
| ScaleRangeUp     | 1.5                   |optional |upper limit of scale factor           |
+------------------+-----------------------+---------+--------------------------------------+
| ShiftRangeLow    | -1.0                  |optional |lower limit of shift factor           |
+------------------+-----------------------+---------+--------------------------------------+
| ShiftRangeUp     |  1.0                  |optional |upper limit of shift factor           |
+------------------+-----------------------+---------+--------------------------------------+
| FixedSyserr      | 1                     |optional |1: fix systematic error; 0: not       |
+------------------+-----------------------+---------+--------------------------------------+
| FixedErrorScale  | 1                     |optional |1: fix error scale; 0: not            |
+------------------+-----------------------+---------+--------------------------------------+
| SyserrRangeLow   | 0.0                   |optional |lower limit of systematic error       |
+------------------+-----------------------+---------+--------------------------------------+
| SyserrRangeUp    | 0.1                   |optional |upper limit of systematic error       |
+------------------+-----------------------+---------+--------------------------------------+
| ErrscaleRangeLow | 0.1                   |optional |lower limit of error scale            |
+------------------+-----------------------+---------+--------------------------------------+
| ErrscaleRangeUp  | 2.0                   |optional |upper limit of error scale            |
+------------------+-----------------------+---------+--------------------------------------+
| SigmaRangeLow    | 1.0e-4                |optional |lower limit of DRW sigma              |
+------------------+-----------------------+---------+--------------------------------------+
| SigmaRangeUp     | 1.0                   |optional |upper limit of DRW sigma              |
+------------------+-----------------------+---------+--------------------------------------+
| TauRangeLow      | 1.0                   |optional |lower limit of DRW tau                |
+------------------+-----------------------+---------+--------------------------------------+
| TauRangeUp       | 1.0e4                 |optional |upper limit of DRW tau                |
+------------------+-----------------------+---------+--------------------------------------+
| FixedCodes       | 1,3                   |optional |codes to be fixed                     |
|                  |                       |         |                                      |
|                  |                       |         |can be multiple codes, use comma (,)  |
|                  |                       |         |to separate, e.g., 1,3                | 
|                  |                       |         |                                      |
|                  |                       |         |this will fix 1st and 3rd code        |
|                  |                       |         |(counting from 0)                     |
+------------------+-----------------------+---------+--------------------------------------+
| FixedScaleCodes  | 1,3                   |optional |codes need to fix the scale           |
|                  |                       |         |                                      |
|                  |                       |         |can be multiple codes, use comma (,)  |
|                  |                       |         |to separate, e.g., 1,3                | 
|                  |                       |         |                                      |
|                  |                       |         |this will fix 1st and 3rd code        |
|                  |                       |         |(counting from 0)                     |
+------------------+-----------------------+---------+--------------------------------------+
| FlagNorm         |  1                    |optional |whether do normalization before       |
|                  |                       |         |intercalibrating                      |
|                  |                       |         |                                      |
|                  |                       |         |1: yes; 0: no                         |
+------------------+-----------------------+---------+--------------------------------------+

Here ``FlagNorm`` specifies whether ``mica2`` does normalization to the light curves of each data codes by their means 
before intercalibrating. This is necessary when there are large offsets between data codes. However, when 
the data has large variability and the time span of data code(s) is short, the means might be seriously biased. 
In this case, doing normalization is not helpful to intercalibration. One may first mannually align the data 
and turn off ``FlagNorm`` option.

After running cali, there is a Python script **plot_for_cali.py** that can used to generate plots,
which generates a PDF file named **PyCALI_results.pdf** and draw a matplotlib 
figure window to show intercalibrated light curves.

Please also refer to :ref:`faq` for more details not covered here.

Python module: pycali
---------------------

An example for using pycali in a Python script is 

.. code-block:: Python
  
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
            fixed_codes=[],
            fixed_scalecodes=[],
            flag_norm=True,
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

Please also refer to :ref:`faq` for more details not covered here.