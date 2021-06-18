function [B] = fun_my_norm(A)
B = sqrt(A(1,:).^2 +A(2,:).^2 + A(3,:).^2);
end

