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
The emssion line fluxes are normalized as

.. math::

  f'_{i,j} = \frac{f_{i, j}}{L_i}\times \frac{L_{i}}{L_0}\frac{C_{i}}{C_0},

where :math:`L i` is the mean of the :math:`i`-th line dataset. This normalization is to enforce that the fluxes are scale with a
same factor as those of continuum. The obtained posterior samples of parameters refer to normalized light curves.
That is to say, one needs to mannually do some convertion to obtain the real parameter values. For scale and shift
parameters,

.. math::

   \varphi \rightarrow \frac{C_0}{C_i} \varphi_i,~~~~~~~~~G_i \rightarrow  C_0 G_i.

For systematic error factor and error scale factors of continuum datasets, 

.. math::

  \epsilon_i \rightarrow C_i \epsilon_i, ~~~~~~~~~~~b_i \rightarrow b_i.

For systematic error factor and error scale factors of line datasets, 

.. math::

  \epsilon_i \rightarrow L'_i \epsilon_i = L_0\frac{C_i}{C_0} \epsilon_i, ~~~~~~~~~~~b_i \rightarrow b_i.