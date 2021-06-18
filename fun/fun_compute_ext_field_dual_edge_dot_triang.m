function Vext=fun_compute_ext_field_dual_edge_dot_triang(NN,nEd,G2,Field_ext,np)

Vext=zeros(nEd,1);

for k=1:nEd
    p_pos=NN(1:3,G2(1,k));
    p_ed=NN(1:3,G2(2,k));
    
    vec1=p_ed-p_pos;
    [PPg_pos,wg1,~]=Gauss_line_nvar2([p_ed,p_pos],np);
    res_1=0;
    for ii = 1:np
        res_1=res_1+0.5*wg1(ii)*dot(Field_ext(PPg_pos(1,ii),PPg_pos(2,ii),PPg_pos(3,ii))',vec1);
    end
    Vext(k)=res_1;
end
    
end


function [PPg,wg,ll]=Gauss_line_nvar2(NN,np)
ll = norm(NN(:,2)-NN(:,1));
[pg,wg]=lgwt(np,-1,1);
pg=pg.';
wg=wg.';
PPg = repmat(NN(:,1),1,np)+repmat(0.5*(NN(:,2)-NN(:,1)),1,np).*repmat((1+pg),3,1);
end