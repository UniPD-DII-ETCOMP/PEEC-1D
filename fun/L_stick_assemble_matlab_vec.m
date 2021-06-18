function [L]=L_stick_assemble_matlab_vec(nSticks,NN,G,radius,npg,Wint,gauss_W,PPghtot,ut_htot,ll_htot)
L=zeros(nSticks,nSticks);
for hh=1:nSticks
    %self-inductance
    PPgh=PPghtot(:,:,hh);
    ll_h=ll_htot(hh);
    aa=log(ll_h/radius(hh)+sqrt((ll_h/radius(hh))^2+1));
    bb=sqrt(1+(radius(hh)/ll_h)^2);
    cc=radius(hh)/ll_h;
    L(hh,hh)=4*ll_h*(aa-bb+cc+(Wint)*0.25);
    ut_h=ut_htot(:,hh);
    %mutual-inductance
    kk=hh+1:nSticks;
        PPk1=NN(1:3,G(1,kk));
        PPk2=NN(1:3,G(2,kk));
        ll_k=ll_htot(kk);
        ut_k=ut_htot(:,kk);
        integ = zeros(length(kk),1);
        for ii=1:npg
            ri=fun_my_norm([PPgh(1,ii)-PPk1(1,:);...
                            PPgh(2,ii)-PPk1(2,:);...
                            PPgh(3,ii)-PPk1(3,:)]);
            rf=fun_my_norm([PPgh(1,ii)-PPk2(1,:);...
                            PPgh(2,ii)-PPk2(2,:);...
                            PPgh(3,ii)-PPk2(3,:)]);
            eps=ll_k./(ri+rf).';
            log_eps=log((1+eps)./(1-eps));
            integ=integ+gauss_W(ii)*log_eps.*fun_my_dot(ut_h,ut_k).';
        end 
        L(hh,kk)=ll_h*integ.';
	    L(kk,hh)=L(hh,kk).';
end 
L=1.0d-7*0.5*L;
end 




