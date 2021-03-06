<p align="center">
	<img src="image2.png" width="850">
</p>

# PEEC 1D 

This directory contains a PEEC code, based on stick elements, for the solution of full-wave electromagnetic problems.
It can be used for research and academic purposes.

Optimized, parallel (OpenMP), fortran90, general versions of this code have been used in, e.g.,

* [R. Torchio et al.,"A 3-D PEEC Formulation Based on the Cell Method for Full-Wave Analyses with Conductive, Dielectric, and Magnetic Media," in IEEE Transactions on Magnetics. doi 10.1109/TMAG.2017.2750319](https://ieeexplore.ieee.org/document/8168416)
* [R. Torchio, "A Volume PEEC Formulation Based on the Cell Method for Electromagnetic Problems From Low to High Frequency," in IEEE Transactions on Antennas and Propagation, vol. 67, no. 12, pp. 7452-7465, Dec. 2019, doi: 10.1109/TAP.2019.2927789](https://ieeexplore.ieee.org/document/8764572)
* [P. Baumgartner et al., "Multi-Objective Optimization of Yagi-Uda Antenna Applying Enhanced Firefly Algorithm With Adaptive Cost Function," in IEEE Transactions on Magnetics. doi: 10.1109/TMAG.2017.2764319](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8168407&isnumber=4479871)
* [T. Bauernfeind et al., "Multi-Objective Synthesis of NFC-Transponder Systems Based on PEEC Method," in IEEE Transactions on Magnetics. doi: 10.1109/TMAG.2017.2771366](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8125565&isnumber=4479871)
* [R. Torchio, P. Bettini and P. Alotto, "PEEC-Based Analysis of Complex Fusion Magnets During Fast Voltage Transients With H-Matrix Compression," in IEEE Transactions on Magnetics. doi: 10.1109/TMAG.2017.2651638](http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7814211&isnumber=7934107)
* [P. Baumgartner et al.,"Synthesis of NFC antenna structure under multi-card condition," in Applied Computational Electromagnetics Society Journal. doi 10.1109/TMAG.2017.2750319](https://ieeexplore.ieee.org/document/8168416)

See the above references for more details and consider citing it.

In this code only conductive media modeled by thin 1D (a.k.a "stick") elements are considered.

-------------------------------------------------------------------

# Description
 
"MAIN_PEEC_1D.m" is the main file you must run to start the code.                      
                                                                                                        
All user-settable quantities, e.g. frequency and resistivity, are contained in the block identified by the 
BEGIN USER SETTINGS / END USER SETTINGS comments.

At the end of the simulation, variable (struct) "sol" contains useful quantities, i.e.:
* sol.I_obj = currents in each stick element of the mesh
* sol.U_obj = electric potentials in each node of the mesh
* sol.Q = electric charges in each node of the mesh
* sol.I_app = currents in each lumped branch
* sol.U_app = electric potentials in each lumped node
* sol.Etot   = total (scattered+external) electric field in target points 
* sol.Esca   = scattered electric field in target points 
* sol.Eext   = external electric field in target points 
-------------------------------------------------------------------

Available test cases
--------------------
Several test cases are contained in separate directories under "test_cases". Each directory contains a description.txt file.
Set the "test_case_dir" variable in "MAIN_PEEC_1D.m" to the appriapriate directory.

User-defined test cases
-----------------------
Follow the instuctions given in "README.txt" inside the "test_cases" directory.

Contacts
-----------------------
Riccardo Torchio (riccardo.torchio@unipd.it)
