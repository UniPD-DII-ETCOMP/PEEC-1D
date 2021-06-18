PEEC 1D (data) by Riccardo Torchio (riccardo.torchio@unipd.it)

How to create a new user-defined test-case:

1. create a new subdirectory
2. the new subdirectory must contain two files:
   a) "description.txt": contains a description of the test case (optional)
   b) "data.mat": contains all data stuctures (see below)
-------------------------------------------------------------------

Description of "data.mat"
(nSticks:number of stick elements, nNodes:number of geometrical points)

"data.mat" MUST contains:
G = 2 x (nSticks:number incidence matrix ( G(1,j)=end node of j-th stick element, G(2,j)=start node of j-th stick element (current flows from start to end))
NN = 3 x nNodes matrix containing coordinates of geometrical points
Lumped = array of struct for lumped circuit elements
    Lumped.node_start=[]; % index of starting node, (if positive it is a node of the mesh, if negative it is a lumped node)
    Lumped.node_end=[];   % index of ending node, (if positive it is a node of the mesh, if negative it is a lumped node)
    Lumped.R=[];          % value of series resistance [Ohm]
    Lumped.L=[];          % value of series self inductance [H]
    Lumped.Cinv=[];       % value of series inverse of self capacitance [F^-1]
    Lumped.value=[];      % value of voltage excitation [V]
Current_Source = array of struct for injected currents 
    Current_Source.node=[]; % index of the node where the current is injected, (if positive it is a node of the mesh, if negative it is a lumped node)
    Current_Source.value=[];% value of the injected current [A]

N.B.: All this variables MUST be contained in "data.mat". For instance,
      even if no Current Sources  are involved in the simulations, the varibale 
      "Current_Source" must be present with empty fields 
      (i.e. Current_Source.node=[];  Current_Source.value=[];)

INFO 
Lumped.node_start
Lumped.node_end
Current_Source.node
can be positive numbers in the range 1:nNodes (nNodes is the number of nodes of the mesh)
when the element is conneted to the device
can be negative numbers if the component is connected to some "appended" circuit node
can be qual to zero if the component is connected to the infinity node

example:
    % FIRST COMPONENT
    Lumped (1).node_start=45; % connected to node 45 (i.e. NN(1:3,45) 
    Lumped (1).node_end=-1;   % connected to the extra appended node 
    Lumped (1).R=50;          % resistance value [Ohm]
    Lumped (1).L=0;           % self inductance value [H]
    Lumped (1).Cinv=0;        % inverse self capacitance value [F^-1]
    Lumped (1).value=3;       % voltage excitation value [V]
    % SECOND COMPONENT (short element)
    Lumped (2).node_start=0;  % connected to the infinity node  
    Lumped (2).node_end=-1;   % connected to the extra appended node
    Lumped (2).R=0;          
    Lumped (2).L=0;
    Lumped (2).Cinv=0;
    Lumped (2).value=0;       % no external excitacion 
% similar for "Current_Source" variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Subdirectories contain some examples of data generation:

"grid": geometrical data are taken from a Comsol mesh file

"other directories": the data are generated in Matlab or uploaded from 
existing .mat files

"line1","line2": can be considered as a easy starting point for
the generation of new user-defined data




 