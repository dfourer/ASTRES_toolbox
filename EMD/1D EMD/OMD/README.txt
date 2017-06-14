PACKAGE EMDOS : EMD by Optimization on Splines, 

Implements the OS method described in [1].



DESCRIPTION
Contains 2 folders:
    - faux/ contains internal function used in the decomposition.
    - figfun/ contains functions which create all the figures of [1]

Contains 2 main functions
    - emdos is the main function for computing the EMD by the OS method. For use make 'help emdos'.
    - paperfig is a script that creates all the figure of [1].


WARNING
This code needs the Spline Toolbox, which is included in the Curve Fitting Toolbox since version 7.11 of Matlab


REFERENCES
  [1] T. Oberlin, S. Meignen and V. Perrier, "An Alternative Formulation
      for the Empirical Mode Decomposition", IEEE Trans. on Signal
      Prcessing, to appear.
  [2] N. Huang, Z. Shen, S. Long, M. Wu, H. Shih, Q. Zheng, N. Yen, C. Tung, and H. Liu, "The empirical mode decomposition
      and the Hilbert spectrum for nonlinear and non-stationary time series analysis," Proceedings of the Royal Society :
      Mathematical, Physical and Engineering Sciences, vol. 454, no. 1971, pp. 903–995, 1998.
  [3] G. Rilling, P. Flandrin, and P. Gonçalvès, "On empirical mode decomposition and its algorithms," in IEEE-EURASIP
      workshop on nonlinear signal and image processing NSIP-03, Grado (I), 2003.




Thomas Oberlin
12.2011
thomas.oberlin@imag.fr
