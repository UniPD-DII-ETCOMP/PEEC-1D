%% ^^^^^^^ Gauss_line_nvar ^^^^^^^^^
function [PPg,ll]=Gauss_line_nvar(NN,gauss_P,np)
ll = norm(NN(1:3,2)-NN(1:3,1));
PPg(1,1:np) = NN(1,1)+0.5*(NN(1,2)-NN(1,1))*(1.0+gauss_P).';
PPg(2,1:np) = NN(2,1)+0.5*(NN(2,2)-NN(2,1))*(1.0+gauss_P).';
PPg(3,1:np) = NN(3,1)+0.5*(NN(3,2)-NN(3,1))*(1.0+gauss_P).';
end