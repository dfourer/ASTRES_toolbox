Prox-EMD

**************************************************************
* author: Nelly Pustelnik                                    *
* institution: Laboratoire de Physique de l'ENS de Lyon      *
* date: Monday, August 25 2014                               *
* License CeCILL-B                                           *
**************************************************************


*****************************************************
* RECOMMENDATIONS:                                  *
* This toolbox is designed to work with             *
* Matlab 7.0                                        *
*****************************************************

------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides a decomposition of signal into a trend and its intrinsic 
mode functions (IMF). The proposed algorithm is an alternative solution 
for the usual Empirical Mode Decomposition by considering optimization tools 
rather than the empirical sifting procedure involved in EMD. 
The proposed method is called Prox-EMD.


This toolbox consists of 4 subfolders:
1) define: MATLAB functions designed for the proposed decomposition algorithm
2) include : contains some standard MATLAB functions
3) EMDs, package_emd, emdos : To compute classical EMD
4) tftb-0.2': Time-frequency toolbox

------------------------------------------------------------------------------------
SPECIFICATIONS for using Prox-EMD:

A demo.m file provides several examples of decomposition:
1) 2 componensts : type        = 'signal_2comp_ex1';
2) 2 componensts : type        = 'signal_2comp_ex2';
3) 2 componensts : type        = 'signal_2comp_ex3';
4) 2 componensts : type        = 'signal_2comp_ex4';
5) 3 componensts : type        = 'signal_3comp_ex1';

The prox_emd.m file is the main function for computing the EMD by the Prox-EMD method. For use make 'help prox_emd'.
An example of a 3-component decomposition (with and without options) is provided in run.m

------------------------------------------------------------------------------------
RELATED PUBLICATIONS:

# N. Pustelnik, P. Borgnat, P. Flandrin 
Empirical Mode Decomposition revisited by multicomponent non smooth convex optimization, 
Signal Processing, Vol. 102, pp. 313--331, Sept. 2014.

# N. Pustelnik, P. Borgnat, and P. Flandrin, 
A multicomponent proximal algorithm for Empirical Mode Decomposition, 
European Signal Processing Conference (EUSIPCO), Bucharest, Romania, August, 27-31, 2012.

http://perso.ens-lyon.fr/nelly.pustelnik/theseE.html
-------------------------------------------------------------------------------------