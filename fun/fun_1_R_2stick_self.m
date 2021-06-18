function [alpha] = fun_1_R_2stick_self(NN,radius)
ll=fun_my_norm(NN(1:3,2)-NN(1:3,1));
aa=log(ll/radius+sqrt((ll/radius)^2+1));
bb=sqrt(1+(radius/ll)^2);
cc=radius/ll;
alpha=2.0*ll*(aa-bb+cc);
end 