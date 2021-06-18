function [A] = fun_A_compute_vec(Nl,Ptar,G1,P0,X_I,normeps,freq)
%%
omega=2*pi*freq;
beta=omega/299792458;
eps0=8.85418781762e-12;
%%
A = zeros(size(Ptar,1),3);
 ii = 1:size(Ptar,1);
    for jj = 1:Nl
        vec_Ri(:,1) = P0(G1(jj,1),1)-Ptar(ii,1);
        vec_Ri(:,2) = P0(G1(jj,1),2)-Ptar(ii,2);
        vec_Ri(:,3) = P0(G1(jj,1),3)-Ptar(ii,3);
        vec_Rf(:,1) = P0(G1(jj,2),1)-Ptar(ii,1);
        vec_Rf(:,2) = P0(G1(jj,2),2)-Ptar(ii,2);
        vec_Rf(:,3) = P0(G1(jj,2),3)-Ptar(ii,3);
        vec_Rm(:,1) = (P0(G1(jj,1),1)+P0(G1(jj,1),1))*0.5-Ptar(ii,1);
        vec_Rm(:,2) = (P0(G1(jj,1),2)+P0(G1(jj,1),2))*0.5-Ptar(ii,2);
        vec_Rm(:,3) = (P0(G1(jj,1),3)+P0(G1(jj,1),3))*0.5-Ptar(ii,3);
        Ri = fun_my_norm(vec_Ri);
        Rf = fun_my_norm(vec_Rf);
        Rm = fun_my_norm(vec_Rm);   
        Ri=Ri+normeps;
        Rf=Rf+normeps;
        ll_k=norm(P0(G1(jj,2),:)-P0(G1(jj,1),:));  
        eps=ll_k./(Ri+Rf);
        log_eps=log((1+eps)./(1-eps));
        ut=P0(G1(jj,2),:)-P0(G1(jj,1),:);   
        ut=ut/ll_k;     
        A(ii,1) = A(ii,1) + 1e-7*X_I(jj)*log_eps.*ut(1).*exp(-1j.*beta.*Rm);
        A(ii,2) = A(ii,2) + 1e-7*X_I(jj)*log_eps.*ut(2).*exp(-1j.*beta.*Rm);
        A(ii,3) = A(ii,3) + 1e-7*X_I(jj)*log_eps.*ut(3).*exp(-1j.*beta.*Rm);
    end    
% end
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



