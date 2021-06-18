clear 
close all
clc

mex -O -largeArrayDims -output L_stick_assemble_for_par L_stick_assemble_par_mex.f90 L_stick_assemble_par.f90
% mex -O -largeArrayDims -output L_edge_vol_assemble_for_par L_edge_vol_assemble_par_mex.f90 L_edge_vol_assemble_par.f90