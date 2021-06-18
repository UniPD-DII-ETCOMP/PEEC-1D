clear 
close all
clc

mex -O -largeArrayDims -output P_stick_assemble_for_par P_stick_assemble_par2_mex.f90 P_stick_assemble_par2.f90
