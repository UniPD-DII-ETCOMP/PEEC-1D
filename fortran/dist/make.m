clear 
close all
clc

% mex -O -largeArrayDims -output dist_edge_edge_for dist_edge_edge_mex.f90 dist_edge_edge.f90
% mex -O -largeArrayDims -output dist_face_face_for dist_face_face_mex.f90 dist_face_face.f90
% mex -O -largeArrayDims -output dist_face_edge_for dist_face_edge_mex.f90 dist_face_edge.f90
mex -O -v -largeArrayDims -output dist_node_node_for dist_node_node_mex.f90 dist_node_node.f90