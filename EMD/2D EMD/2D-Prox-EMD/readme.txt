Toolbox_PHT_2D

*******************************************************
* author: Jérémy Schmitt                              *
* institution: Laboratoire de Physique de l'ENSL      *
* date: Wednesday, October 29 2014     	              *
* License CeCILL-B                                    *
*******************************************************


*****************************************************
* RECOMMENDATIONS:                                  *
* This toolbox is designed to work with              *
* Matlab 7.0                                        *
*****************************************************

------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides a spectral analysis of a 2D image called Prony-Huang transform. This method is based on the principle of Hilbert-Huang transform and involves variational methods.
The procedure is based on two steps:
- A mode decomposition step decomposing the image into a trend and its intrinsic mode functions (IMF). The proposed algorithm is an alternative solution for the usual 2D Empirical Mode Decomposition and involves optimization tools rather than the usual sifting process involved in EMD. This method is called Variational 2D-EMD.
- A spectral analysis step performed patch-wise on the IMFs, based on 2D Prony method.


This toolbox consists of 2 subfolders:
1) include: MATLAB functions designed for the proposed algorithm
2) colormap_utilities: utilities for subplots in hsv colormaps (useful for the plotting of Prony results)

------------------------------------------------------------------------------------
SPECIFICATIONS for using Prox-EMD:

Two demos files are proposed. demo_emd.m provides examples of Variational 2D-EMD (without the spectral analysis step). demo_pht.m provides examples of Prony-Huang spectral analysis.

Both demo files provide several examples of decomposition:
1) `example1’: 2 sinusoïdal components
2) `example2’: 3 sinusoïdal components
3) ‘example3’: simulations from IEEE-TIP paper : 2 FM components and 2 piecewise constant patches

The main functions are run_prox_emd_2D.m and run_prony_2D.m.
run_prox_emd_2D.m performs Variational 2D-EMD on an image. run_prony_2D.m performs Prony spectral analysis on each IMF obtained with run_prox_emd_2D.

------------------------------------------------------------------------------------
RELATED PUBLICATIONS:

# J. Schmitt, N. Pustelnik, P. Borgnat, P. Flandrin, L. Condat
2-D Prony–Huang Transform: A New Tool for 2-D Spectral AnalysisIEEE Transactions on Image Processing, Vol. 23, Issue 12, pp. 5233—5248, Dec. 2014

# J. Schmitt, N. Pustelnik, P. Borgnat, P. Flandrin
2-D Hilbert-Huang TransformInternational Conference on Acoustics, Speech and Signal Processing (ICASSP), Florence, Italy, May, 4-9, 2014
-------------------------------------------------------------------------------------