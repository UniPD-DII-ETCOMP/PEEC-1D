function [A_eapp,A_napp,R_app,L_app,Cinv_app,Us_app]=set_Lumped(nn_dual,n_Vsource,nn_app,Voltage_Source)

A_eapp = zeros(nn_dual,n_Vsource);
A_napp = zeros(nn_app,n_Vsource);
R_app = zeros(n_Vsource,n_Vsource);
L_app = zeros(n_Vsource,n_Vsource);
Cinv_app = zeros(n_Vsource,n_Vsource);
Us_app = zeros(n_Vsource,1);

for k=1:n_Vsource
    %Lumped element connected from INF node to OBJ node
    if sign(Voltage_Source(k).node_end) == 1 && sign(Voltage_Source(k).node_start) == 0
        A_eapp(Voltage_Source(k).node_end,k)=1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;
    %Lumped element connected form INF node to APPENDED node    
    elseif sign(Voltage_Source(k).node_end) == -1 && sign(Voltage_Source(k).node_start) == 0
        A_napp(abs(Voltage_Source(k).node_end),k)=1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;
    %Lumped element connected from APPENDED node to OBJ node
    elseif sign(Voltage_Source(k).node_end) == 1 && sign(Voltage_Source(k).node_start) == -1
        A_eapp(Voltage_Source(k).node_end,k)=1;
        A_napp(abs(Voltage_Source(k).node_start),k)=-1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;
    %Lumped element connected from OBJ node to APPENDED node
    elseif sign(Voltage_Source(k).node_end) == -1 && sign(Voltage_Source(k).node_start) == 1
        A_eapp(Voltage_Source(k).node_start,k)=-1;
        A_napp(abs(Voltage_Source(k).node_end),k)=1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;
    %Lumped element connected from OBJ node to OBJ node
    elseif sign(Voltage_Source(k).node_end) == 1 && sign(Voltage_Source(k).node_start) == 1
        A_eapp(Voltage_Source(k).node_start,k)=-1;
        A_eapp(Voltage_Source(k).node_end,k)=1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;
    %Lumped element connected from APPENDED node to APPENDED node
    elseif sign(Voltage_Source(k).node_end) == -1 && sign(Voltage_Source(k).node_start) == -1
        A_napp(abs(Voltage_Source(k).node_start),k)=-1;
        A_napp(abs(Voltage_Source(k).node_end),k)=1;
        R_app(k,k)=Voltage_Source(k).R;
        L_app(k,k)=Voltage_Source(k).L;
        Cinv_app(k,k)=Voltage_Source(k).Cinv;
        Us_app(k)=Voltage_Source(k).value;        
    end  		
end  

A_eapp = sparse(A_eapp);
A_napp = sparse(A_napp);

end