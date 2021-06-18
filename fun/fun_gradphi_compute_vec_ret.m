function [E] = fun_gradphi_compute_vec_ret(Ptar,P0,q,freq)
Nl=length(q);
%%
omega=2*pi*freq;
beta=omega/299792458;
eps0=8.85418781762e-12;
%%
E = zeros(size(Ptar,1),3);
 ii = 1:size(Ptar,1);
    for jj = 1:Nl
        vec_R(:,1) = -P0(jj,1)+Ptar(ii,1);
        vec_R(:,2) = -P0(jj,2)+Ptar(ii,2);
        vec_R(:,3) = -P0(jj,3)+Ptar(ii,3);
        R = fun_my_norm(vec_R);         
        E(ii,1) = E(ii,1) + (1/(4*pi*eps0)).*q(jj)*(exp(-1j.*beta.*R)).*(...
                             vec_R(:,1)./R.^3+1j*omega.*vec_R(:,1)./(299792458*R.^2));
        E(ii,2) = E(ii,2) + (1/(4*pi*eps0)).*q(jj)*(exp(-1j.*beta.*R)).*(...
                             vec_R(:,2)./R.^3+1j*omega.*vec_R(:,2)./(299792458*R.^2));
        E(ii,3) = E(ii,3) + (1/(4*pi*eps0)).*q(jj)*(exp(-1j.*beta.*R)).*(...
                             vec_R(:,3)./R.^3+1j*omega.*vec_R(:,3)./(299792458*R.^2));          
    end    
end
function [B] = fun_my_norm(A)
B = sqrt(A(:,1).^2 +A(:,2).^2 + A(:,3).^2);
end
function [C] = fun_my_dot(A,B)
C = A(1)*B(1)+A(2)*B(2)+A(3)*B(3);
end
function [C] = fun_my_cross(A,B)
C = [A(2)*B(3)-A(3)*B(2),A(3)*B(1)-A(1)*B(3),A(1)*B(2)-A(2)*B(1)];
end



