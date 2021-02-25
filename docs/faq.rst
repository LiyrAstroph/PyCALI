.. _faq:

**************************
Frequently Asked Questions
**************************

1. **The fluxes are all in physical units not magnitudes, right?**
  
   Yes.

2. **Does any of the light curves need to overlap with the other light curves?**
   
   In principle the code does not require that.

   However, it would be better that the different datasets have some 
   overlap periods, otherwise, the code just effectively align 
   the means of all datasets with the reference dataset.   
   It is worth mentioning that different datasets DO NOT need 
   to have exactly same epochs.

3. **Should the codes of continuum and emission-line light curves be same?**
   
   Yes. The orders of the codes appearing in the input data also should be exactly same.

   A dataset is permitted to have zero points.

4. **Can PyCALI only merge continuum light curves?**
   
   Yes. Switch off the option for lines.

5. **Can PyCALI only merge line light curves?**

   Yes. One can treat line light curves as continuum and fix the shift factor (to zero).  