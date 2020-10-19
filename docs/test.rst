*******
Tests
*******

A mock light curve is generated using the damped random walk model and then passed to pyCALI. 
The obtained estimates for scale factors, shift factors, and systematic errors are consistent 
with the input values.

.. figure:: _static/test_lightcurve.jpg
  :scale: 30 %
  :align: center
  
  A mock light curve with five datasets and the intercalibration.

.. figure:: _static/test_syserr.jpg
  :scale: 30 %
  :align: center
  
  The obtained posterior distributions of the systematic error factors. Blue lines represent the input values.
