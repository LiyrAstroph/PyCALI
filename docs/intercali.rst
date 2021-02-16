************************************************
Intercalibration: Methodology and Implementation 
************************************************

Basic Equations
===============
By taking :math:`5100` A continuum flux densities and H :math:`\beta` emission line fluxes 
as an example, the intercalibration between datasets observed 
by different telescopes is written

.. math::
  
  F_\lambda (5100~\text{A}) = \varphi \cdot F_\lambda (5100~\text{A})_{obs} - G,

  F(H\beta) = \varphi \cdot F(H\beta)_{obs},

where :math:`\varphi` is a multiplicative factor and :math:`G` is an additive factor.

Reform the above equations with vectors. We now have :math:`m` measurements of a true signal :math:`y` by :math:`k` different telescopes

.. math::
  
  {f}_{obs} = {\Phi}^{-1}({y} + {G}) + {bn} + {\epsilon},

  {f}_{cali} = {\Phi}({f}_{obs}+ {bn} + {\epsilon}) - {G},

where :math:`n` are reported measurement noises, :math:`\epsilon` are unknown systematic errors, 
:math:`b` is an :math:`m\times m` diagnoal matrix for error scale. The above equations stands for 
continuum flux densities. For emission line fluxes, :math:`G=0`.

PyCALI uses the damped random walk process to describe the variations and employs a Bayesian 
framework to determine the best estimates for the parameters by exploring the posterior probability distribution
with the diffusive nested sampling algorithm.

Implementation
==============
In real implementation, the first dataset is set to be the reference with :math:`\varphi=1` and :math:`G=0`. 
The input light curve of each dataset is first normalized before being passed to intercalibration. The continuum fluxes are normalized as  

.. math::
  
  f'_{i,j} = \frac{f_{i, j}}{C_i},

where :math:`f_{i, j}` is the :math:`j`-th point of the :math:`i`-th continuum dataset and :math:`C_i` is the mean. 
The emssion line fluxes are normalized as

.. math::

  f'_{i,j} = \frac{f_{i, j}}{L_i}\times \frac{L_{i}}{L_0}\frac{C_{i}}{C_0} =  \frac{f_{i, j}}{L_0\frac{C_{i}}{C_0}} = \frac{f_{i, j}}{L'_i},

where :math:`L_i` is the mean of the :math:`i`-th line dataset, and 

.. math::

  L'_i = L_0\frac{C_{i}}{C_0}.

This normalization is to enforce that the fluxes are scaled with a
same factor as those of continuum. The obtained posterior samples of parameters refer to normalized light curves.
That is to say, one needs to mannually do some convertions to obtain the real parameter values. For scale and shift
parameters,

.. math::

   \varphi_i \rightarrow \frac{C_0}{C_i} \varphi_i,~~~~~~~~~G_i \rightarrow  C_0 G_i.

For systematic error factor and error scale factors of continuum datasets, 

.. math::

  \epsilon_i \rightarrow C_i \epsilon_i, ~~~~~~~~~~~b_i \rightarrow b_i.

For systematic error factor and error scale factors of line datasets, 

.. math::

  \epsilon_i \rightarrow L'_i \epsilon_i = L_0\frac{C_i}{C_0} \epsilon_i, ~~~~~~~~~~~b_i \rightarrow b_i.

Outputs
=======
PyCALI generates a number of outputs to the directory **./data/**. If this directory does not exist, PyCALI will create it automatically.
Main output files are 

* xxx_cali

  intercalibrated light curve, **xxx** represents the name of the input data.

* posterior_sample.txt

  posterior samples for intercalibration parameters. The columns are: 
  sigmad and taud (damped random walk model parameters), scale factors, shift factors, systematic error factors, 
  and error scale factors. (Note that the values of these parameters refer to the light curves normalized by their means, see above).

* factor.txt 

  The estimated values of intercalibration parameters, determined from the means and standard deviations of the posterior samples.
  One may do more sophiciated statitics using the posterior samples.
  (Note again that these values refer to the light curves normalized by their means, see above).

* xxx_recon.txt
  
  reconstructions to the intercalibrated light curves using the dampled random walk model.

Special Notes
=============

* When a dataset has number of continuum points less than or equal to 2 and no line points, the shift factor is fixed to zero.

* The scale and shift factors are highly degenerated. PyCALI implicitly take this degenecy into account when 
  peforming Bayesian sampling.

* PyCALI employs the diffusive nested sampling package CDNest (https://github.com/LiyrAstroph/CDNest) to generate 
  posterior samples of intercalibration parameters. The diffusive nested sampling algorithm is developed by Brewer et al. (2011; 
  https://github.com/eggplantbren/DNest4).

Future Improvements
===================

Presently, PyCALI relies on the damped random walk process to describe light curves. An adaption to flexible variability models 
is highly worthwhile.