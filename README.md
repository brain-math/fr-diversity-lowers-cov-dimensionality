# fr-diversity-lowers-cov-dimensionality

Code accompanying the paper 
"_Neuronal firing rate diversity lowers the dimension of population covariability_" (2026) by Gengshuo John Tian, Ou Zhu, Vinay Shirhatti, Charles M. Greenspon, John E. Downey, David J. Freedman, and Brent Doiron.

`Sim.m` contains code for numerical simulations (Figs. 2, 3a, 6c, S3, S8). 
Each section can be run independently to obtain the figures it is responsible for.
Code for parameter sweeps are omitted since they take too much time and can be obtained by straight-forward modifications of the code for single parameters if desired.

`Analysis.m` illustrates the data analysis pipeline (Fig. 3d and Fig. 4). 
Requires public V1-V2 data from Zandvakili & Kohn (2015), which can be downloaded via https://crcns.org/data-sets/vc/v1v2-1/about_v1v2-1

The core functions are:
- `computepdf.m` computes the eigenvalue distribution by solving for the Cauchy transform using Eq. (S9).
- `prdat.m` computes the temporal sampling-corrected participation ratio D(C) using Eq. (6).
- `fitkng.m` fits the best model for the data eigenvalue distribution, which is then used by `detectol.m` to detect outlier eigenvalues.
