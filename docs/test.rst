*******
Tests
*******

A mock light curve is generated using the damped random walk model and then can be passed to PyCALI. 
PyCALI provides a function to generate mock light curves

.. code-block:: python

  pycali.generate_mock_data()

This will generate mock continnum and emission-line light curves (**sim_cont.txt** and **sim_line.txt**) 
and place them to the directory **./data**. An example Python script **example_mock.py** in the source code 
are provided to do tests with mock data. 

The obtained estimates for scale factors, shift factors, and systematic errors are shown below, which 
are generally consistent with the input values.

.. figure:: _static/test_lightcurve.jpg
  :scale: 25 %
  :align: center
  
  A mock light curve with five datasets and the intercalibration.

.. figure:: _static/test_syserr.jpg
  :scale: 25 %
  :align: center
  
  The obtained posterior distributions of the systematic error factors. Blue lines represent the input values.
