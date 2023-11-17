Companion page for Graphic Equalizer Design with Symmetric Biquad Filters

http://research.spa.aalto.fi/publications/papers/waspaa19-sge/

The MATLAB files in this web folder are related to the following conference article:

J. Liski, J. Ramo, and V. Valimaki, "Graphic equalizer design with symmetric biquad 
filters," in Proceedings of IEEE Workshop on Applications of Signal Processing to Audio 
and Acoustics (WASPAA), New Paltz, NY, USA, 2019.

MATLAB code for Figures 1, 2, 3, and 4 of the above paper are included (Fig1.m, Figs2.m, 
Figs3.m, Figs4.m). These scripts can be executed separately, but please note that they use
some of the additional MATLAB files, which must be made available for MATLAB.

Additional MATLAB files:
- acge_noIter.m = implements an older state-of-the-art graphical EQ without iteration step
- ho_geq_SGE.m = implements a higher-order cascade EQ utilizing fourth-order sections
- interactionMatrix2.m = function used to calculate the interaction matrix for ACGE
- interactionMatrixSGE.m = function used to calculate the interaction matrix for SGE
- pareq.m = function that creates parametric EQ filters according to S. Orfanidis, 
Introduction to Signal Processing, p. 594.
- peq.m = function that creates parametric EQ filters with adjustable  Nyquist gain 
according to S. Orfanidis, "Digital parametric equalizer design with prescribed Nyquist-
frequency gain"
- peq_SGE.m = proposed modification to the peq.m function
- sge.m = design graphic EQ with symmetric biquad filters

In Espoo, Finland, 17 October 2019

Juho Liski

Aalto University, Dept. of Signal Processing and Acoustics

(end of file)