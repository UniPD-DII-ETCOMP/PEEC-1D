function [log_eps] = fun_1_R_stick_vec(PP1,PP2,PP0,ll)
ri=fun_my_norm([PP0(1)-PP1(1,:);...
                PP0(2)-PP1(2,:);...
                PP0(3)-PP1(3,:)]);
rf=fun_my_norm([PP0(1)-PP2(1,:);...
                PP0(2)-PP2(2,:);...
                PP0(3)-PP2(3,:)]);
%ll=fun_my_norm(PP(1:3,2)-PP(1:3,1))
eps=ll./(ri+rf);
log_eps=log((1+eps)./(1-eps));
end 
%%
% fun_1_R_stick