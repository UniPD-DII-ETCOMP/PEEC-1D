function [P] = P_stick_assemble_matlab_vec(nNodes,NN,G,radius,npg,ne_cap_max,Cap_Elem,...
               NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,gauss_W)
%
%
epsilon0=8.85418781762d-12;
%%
P=zeros(nNodes,nNodes);
%%
for hh=1:nNodes
    ne_cap_hh=Cap_Elem(1,hh);
    idE_hh=Cap_Elem(2:ne_cap_max+1,hh);
    P_self_jj_jj=0.0d0;
    P_self_jj_ii=0.0d0;
    ll_tot_hh=0.0d0;
    for jj=1:ne_cap_hh
        %self-local-capacitance
        NN_edge=NN(1:3,G(1:2,idE_hh(jj)));
        NN_jj(1:3,1) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)); 
        NN_jj(1:3,2) = NN(1:3,hh);
        alpha_jj= fun_1_R_2stick_self(NN_jj,radius(idE_hh(jj)));
        P_self_jj_jj=P_self_jj_jj+alpha_jj; 
%        call Gauss_line_nvar(NN_jj,gauss_P,gauss_W,npg,PPg_jj,whg_jj,ll_jj)
%         [PPg_jj,ll_jj] = Gauss_line_nvar(NN_jj,gauss_P,npg);
        PPg_jj=PPg_xx(1:3,1:npg,jj,hh);
        ll_jj=ll_xx(jj,hh); 
        ll_tot_hh=ll_tot_hh+ll_jj;
		%mutual-local-capacitance
        for ii=jj+1:ne_cap_hh
            NN_ii(1:3,1) = reshape(NN_xx1(1:3,ii,hh),3,1); 
            NN_ii(1:3,2) = NN(1:3,hh);
            ll_ii=ll_xx(ii,hh);
            %gauss
            integ_self=0.0d0;
            for ff=1:npg
                log_eps = fun_1_R_stick(NN_ii,PPg_jj(1:3,ff),ll_ii);
                integ_self=integ_self+gauss_W(ff)*log_eps;
            end 
            P_self_jj_ii=P_self_jj_ii+0.5*ll_jj*integ_self;
        end 
    end 
	P(hh,hh)=(P_self_jj_jj+2*P_self_jj_ii)/(4*pi*epsilon0*ll_tot_hh^2);
end
% return
for hh=1:nNodes
    ne_cap_hh=Cap_Elem(1,hh); 
    ll_tot_hh=ll_tot_xx(hh).';
    kk=hh+1:nNodes;
    
       len_zz=length(kk);
	   P_mutual_jj_ii = zeros(1,len_zz);
       ll_tot_kk=ll_tot_xx(kk).';
       
	   for jj = 1:ne_cap_hh
          PPg_jj=PPg_xx(1:3,1:npg,jj,hh);
          ll_jj=ll_xx(jj,hh); 
          for ii = 1:ne_cap_max
             NN_ii1 = reshape(NN_xx1(1:3,ii,kk),3,len_zz); 
             NN_ii2 = reshape(NN_xx2(1:3,ii,kk),3,len_zz);%NN(1:3,kk);
             ll_ii=ll_xx(ii,kk);		 
			 integ_mutual=0.0;
			 for ff=1:npg
                 [log_eps] = fun_1_R_stick_vec(NN_ii1,NN_ii2,PPg_jj(1:3,ff),ll_ii);
                  integ_mutual=integ_mutual+gauss_W(ff).*log_eps;			 
			 end 
	         P_mutual_jj_ii=P_mutual_jj_ii+0.5*ll_jj.*integ_mutual;
		  end  
		  %glob_P=glob_P+P_mutual_jj_ii
	   end   
	   P(hh,kk)=P_mutual_jj_ii./(ll_tot_kk*ll_tot_hh*4*pi*epsilon0);
	   P(kk,hh)=P(hh,kk).';
%     end
end
%%
end
%


