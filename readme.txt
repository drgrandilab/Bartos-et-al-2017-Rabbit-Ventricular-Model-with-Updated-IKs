Matlab code of the Bartos et al. model of IKs in rabbit ventricular myocytes

___________________________________________________________________________________________
Contents:

readme.txt				  this file

rabbit_myofilament_masterCompute.m	loads initial conditions and runs the simulation
rabbit_myofilament_masterODEfile.m	integrates the following model components
rabbit_myofilament_barODEfile.m		  beta-adrenergic (PKA) phosphorylation module
rabbit_myofilament_camkiiODEfile.m	CaMKII phosphorylation module
rabbit_myofilament_camODEfile.m		  CaM module
rabbit_myofilament_eccODEfile.m		  excitation-contraction coupling module

AP_Bartos.mat				voltage signal used in AP-clamp simulations

Other .mat files			      initial conditions
- yfin_iks_isometric_1Hz		steady-state @ 1 Hz, no ISO
- yfin_iks_isometric_3Hz		steady-state @ 3 Hz, no ISO
- yfin_iks_isometric_3Hz_0p020iso_240s	3 Hz, 240 s ISO administration (20 nM)
- yfin_iks_isometric_3Hz_0p050iso_240s	3 Hz, 240 s ISO administration (50 nM)
___________________________________________________________________________________________


This model has been obtained by updating the formulation of IKs present in the excitation-
contraction coupling module in the parent model developed by Negroni et al.


References:

Bartos DC, Morotti S, Ginsburg KS, Grandi E, Bers DM.
Quantitative analysis of the Ca2+-dependent regulation of delayed rectifier K+ current IKs
in rabbit ventricular myocytes.
J Physiol. 2017 Apr 1;595(7):2253-2268. doi: https://doi.org/10.1113/JP273676

Negroni JA, Morotti S, Lascano EC, Gomes AV, Grandi E, Puglisi JL, Bers DM.
ß-adrenergic effects on cardiac myofilaments and contraction in an integrated rabbit
ventricular myocyte model.
J Mol Cell Cardiol. 2015 Apr;81:162-75. doi: https://doi.org/10.1016/j.yjmcc.2015.02.014

Please cite the above papers when using this model.
