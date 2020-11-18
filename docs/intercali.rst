************************************************
Intercalibration: Methodology and Implementation 
************************************************

Basic Equations
===============
By taking :math:`5100\AA` continuum flux densities as examples, the intercalibration between datasets observed 
by different telescopes is written

.. math::
  
  F_\lambda (5100~\text{\AA}) = \varphi \cdot F_\lambda (5100~\text{\AA})_{obs} - G,

where :math:`\varphi` is a multiplicative factor and :math:`G` is an additive factor.

Reform the above equations with vectors. We now have :math:`m` measurements of a true signal :math:`y` by :math:`k` different telescopes

.. math::
  
  {f}_{obs} = {\Phi}^{-1}({y} + {G}) + {bn} + {\epsilon},

  {f}_{cali} = {\Phi}({f}_{obs}+ {bn} + {\epsilon}) - {G},

where :math:`n` are reported measurement noises, :math:`\epsilon` are unknown systematic errors, 
:math:`b` is an :math:`m\times m` diagnoal matrix for error scale.

pyCALI uses the damped random walk process to describe the variations and employs a Bayesian 
framework to determine the best estimates for the parameters by exploring the posterior probability distribution
with the diffusive nested sampling algorithm.

Implementation
==============
In real implementation, the first dataset is set to be the reference with :math:`\varphi=1` and :math:`G=0`. 
The input light curve of each dataset is first normalized before being passed to intercalibration. The continuum fluxes are normalized as  

.. math::
  
  f'_{i,j} = \frac{f_{i, j}}{C_i},

where :math:`f_{i, j}` is the :math:`j`-th point of the :math:`i`-th continuum dataset and :math:`C_i` is the mean. 