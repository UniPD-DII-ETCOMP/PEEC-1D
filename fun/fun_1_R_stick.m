function [log_eps] = fun_1_R_stick(PP,PP0,ll)
ri=fun_my_norm(PP0-PP(1:3,1));
rf=fun_my_norm(PP0-PP(1:3,2));
%ll=fun_my_norm(PP(1:3,2)-PP(1:3,1))
eps=ll/(ri+rf);
log_eps=log((1+eps)/(1-eps));
end 