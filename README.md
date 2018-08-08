# HMC
All files are writted in MATLAB code.

1) The file synthetic_data_generate.m generates synthetic data with inpus "d" (the dimension) and "r" (the number of data points).

2) The file HMC_and_Langevin.m  runs the HMC or Langevin algorithms.  Which algorithm is run (HMC or Langevin) and the parameters "eta" and "T" of the algorithms can be set by the user by changing these parameters in the .m file.


3) The file L2__Vs__L_infty.m  approximates the (rescaled) values of the Lipschitz Hessian constant $L2$ and the infinity-norm Lipschitz constant $L_\infty$ at different values of the dimension $d$ for synthetic data.
