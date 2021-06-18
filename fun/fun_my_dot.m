%%
function [C] = fun_my_dot(A,B)
C = A(1)*B(1,:)+A(2)*B(2,:)+A(3)*B(3,:);
end