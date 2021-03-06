clc
clear
close all
restoredefaultpath
%% BEGIN USER SETTINGS
% problem definition
test_case_dir = 'nfc1';
plotflag = 1; %graphics postprocessing flag (1 = yes, 0 = no)
freq = 3e8; % set the value of the frequency [Hz]
radius = 0.1e-3; % set the radius of the wires  [m]
rho = 1/57e6;% set the value of the resistivity [Ohm m]
% algorithmic settings
ret = 1; %consider retarded potentials (1 = yes, 0 = no)
Wint = 0; %consider internal energy of wires for the self inductance computation (1 = yes, 0 = no)
np_gauss = 2; % number of gauss points for the numerical integration
E_external = 0; % external electric field (1 = yes, 0 = no)
% External field (only active if E_external==1)
E_ext=@(x,y,z) [(y-0.5),-(x-0.5),0.*z]; % set the external electric field distribution
% number of thread
N.thread=22;
% E field computation in target points 
E_target  = 1; % extrnal field E is computed in target points
TP =[0 0 0;
     0 0 1;
     0 0 2;
     0 0 3;
     0 0 4;
     0 0 5;
     0 0 6;
     0 0 7;
     0 0 8;
     0 0 9]*1e-3;%*1e-3;  
% END USER SETTINGS
%% data input (load data.mat examples from "test_case_dir" directory in "data_generation" directory)
cd test_cases
cd(test_case_dir)
load data.mat
try 
type description.txt;
catch
    
end
disp(' ')
cd ..
cd ..
%% computing matrices
disp('Pre-processing ...')
tic
Matrix_P0=NN.';
xmin=min(NN(1,:))-1e-3;ymin=min(NN(2,:))-1e-3;zmin=min(NN(3,:))-1e-3;xmax=max(NN(1,:))+1e-3;ymax=max(NN(2,:))+1e-3;zmax=max(NN(3,:))+1e-3;
bar_stick = 0.5*(Matrix_P0(G(1,:),:)+Matrix_P0(G(2,:),:));
eps0=8.854e-12;
mu0=4*pi*1e-7;
cd('fun'); addpath(genpath(pwd)); cd ..
cd('fortran'); addpath(genpath(pwd)); cd ..
omega=freq*2*pi;
omega2=omega;
if omega == 0
    omega2=1;
    warning('lumped capacitance are ignored')
end
beta = omega/(1/sqrt(eps0*mu0));
nSticks = size(G,2);
nNodes = size(NN,2);
lgt = compute_length(NN,G);
val_neg_g=G(1,:);
val_pos_g=G(2,:);
Aobj=sparse([1:nSticks],val_neg_g,-ones(nSticks,1),nSticks,nNodes);
Aobj=Aobj+sparse([1:nSticks],val_pos_g,ones(nSticks,1),nSticks,nNodes);
Aobj=Aobj.';
% Appended elements 
N.lumped=size([Lumped.value],2);
N.cur_sou=size([Current_Source.node],2);
app=[[Lumped.node_start],[Lumped.node_end],[Current_Source.node]];
app=unique(app);
app=find(app<0);
N.node_app=length(app);
toc
disp(' ')
disp('INFO:')
disp(['   ',num2str(nSticks) ' sticks'])
disp(['   ',num2str(nNodes) ' nodes'])
disp(['   ',num2str(N.lumped) ' lumped edges'])
disp(['   ',num2str(N.node_app) ' lumped nodes'])
disp(['   ',num2str(N.cur_sou) ' injected currents'])
disp(' ')
%%
disp('Computing matrices ...')
tic
% R 
R = lgt.*(rho./(pi*radius.^2));
R = diag(R);
R = sparse(R);
% P 
ne_cap_max = max(full(sum(abs(Aobj),2)));
Cap_Elem = zeros(1+ne_cap_max,nNodes);
for k=1:nNodes
    idE = find(Aobj(k,:)); %edges connected to node k
    nE = length(idE);
    Cap_Elem(1,k) = nE;
    Cap_Elem(2:nE+1,k) = idE;
end
ne_cap_tot = sum(Cap_Elem(1,:));
try % fortran
    P = P_stick_assemble_for_par(nSticks,nNodes,NN,G,radius*ones(1,nSticks),np_gauss,ne_cap_max,Cap_Elem,N.thread);
    if ret == 1
        dist_P = dist_node_node_for(NN,nNodes,N.thread);
        P = P.*exp(-1j*dist_P*beta);
    end
catch % matlab
    warning('- MEX function P not supported, try to re-mex it: run /fortran/make.m. Switched to slow Matlab function.')
    [gauss_P,gauss_W]=lgwt(np_gauss,-1,1);
    NN_xx1=zeros(3,ne_cap_max,nNodes);
    NN_xx2=zeros(3,ne_cap_max,nNodes);
    PPg_xx=zeros(3,np_gauss,ne_cap_max,nNodes);
    ll_xx=zeros(ne_cap_max,nNodes);
    ll_tot_xx=zeros(nNodes,1);
    for hh=1:nNodes 
        ne_cap_hh=Cap_Elem(1,hh); % number of elements 
        idE_hh=Cap_Elem(2:ne_cap_max+1,hh);	
        glob_P=0.0;
       for jj = 1:ne_cap_hh
          NN_edge=NN(1:3,G(1:2,idE_hh(jj))); % take end points
          NN_xx1(1:3,jj,hh) = 0.5*(NN_edge(1:3,1)+NN_edge(1:3,2)); % build capacitive element
          NN_xx2(1:3,jj,hh) = NN(1:3,hh); 
          [PPg_xx(1:3,1:np_gauss,jj,hh),ll_xx(jj,hh)] = Gauss_line_nvar([NN_xx1(:,jj,hh),NN_xx2(:,jj,hh)],gauss_P,np_gauss); %gauss points
          ll_tot_xx(hh)=ll_tot_xx(hh)+ll_xx(jj,hh); % length update
       end
    end    
    P=P_stick_assemble_matlab_vec(nNodes,NN,G,radius*ones(1,nSticks),np_gauss,ne_cap_max,Cap_Elem,...
                 NN_xx1,NN_xx2,PPg_xx,ll_xx,ll_tot_xx,gauss_W);
    if ret == 1
        dist_P = dist_matlab(NN);
        P = P.*exp(-1j*dist_P*beta);
    end
end
% L
try % fortran
    L = L_stick_assemble_for_par(nSticks,nNodes,NN,G,radius*ones(1,nSticks),np_gauss,Wint,N.thread);
    if ret == 1
        dist_L = dist_edge_edge_for(NN,G,nSticks,nNodes,N.thread);
        L = L.*exp(-1j*dist_L*beta);
    end
catch % matlab
    warning('- MEX function L not supported, try to re-mex it: run /fortran/make.m. Switched to slow Matlab function.')
    [gauss_P,gauss_W]=lgwt(np_gauss,-1,1); % gauss points and weights (local)
    PPgh=zeros(3,np_gauss,nSticks);
    ll_h=zeros(nSticks,1);
    ut_h=zeros(3,nSticks);
    for ii = 1:nSticks
        PPh=NN(1:3,G(1:2,ii)); 
        [PPgh(:,:,ii),ll_h(ii)] = Gauss_line_nvar(PPh,gauss_P,np_gauss); % gauss points global
        ut_h(:,ii)=(PPh(1:3,2)-PPh(1:3,1))/ll_h(ii);
    end   
    L=L_stick_assemble_matlab_vec(nSticks,NN,G,radius*ones(1,nSticks),np_gauss,Wint,gauss_W,PPgh,ut_h,ll_h);
    if ret == 1
        dist_L = dist_matlab(bar_stick.');
        L = L.*exp(-1j*dist_L*beta);
    end
end
% rhs
if E_external==1
Us = fun_compute_ext_field_dual_edge_dot_triang(NN,nSticks,G,E_ext,7);
else
Us = zeros(nSticks,1);   
end
Is = zeros(nNodes,1);
Is_napp = zeros(N.node_app,1);
% app
[A_eapp,A_napp,R_app,L_app,Cinv_app,Us_app] = set_Lumped(nNodes,N.lumped,N.node_app,Lumped); 
for k=1:N.cur_sou 
    if sign(Current_Source(k).node) == 1  %current generator connected to a OBJ node
        Is(Current_Source(k).node)=Current_Source(k).value;
    elseif sign(Current_Source(k).node) == -1  %current generator connected to a APPENDED node
        Is_napp(abs(Current_Source(k).node))=Current_Source(k).value;
    end
end
toc
disp(' ')
%%
disp('Build system M ...')
tic
M = [ R+1j*omega*L, ... (e_obj/e_obj)
      Aobj.',... % (e_obj/n_obj)
      sparse(nSticks,N.lumped), ... (e_obj/e_app)
      sparse(nSticks,N.node_app); ... (e_app/n_app)  
      P*Aobj, ... (n_obj/e_obj)
      -(1j*omega)*speye(nNodes,nNodes),...  (n_obj/n_obj)
      P*A_eapp, ... (n_obj/e_app)
      sparse(nNodes,N.node_app); ... (n_obj/n_app)
      sparse(N.lumped,nSticks), ... (e_app/e_obj)
      A_eapp.', ... (e_app/n_obj)
      R_app+1j*omega*L_app+(1/(1j*omega2))*Cinv_app, ...(e_app/e_app)   (Cinv deleted for DC simulations)
      A_napp.'; ... (e_app/n_app)
      sparse(N.node_app,nSticks), ... (n_app/e_obj)
      sparse(N.node_app,nNodes),... (n_app/v_obj)
      (1j*omega2)*A_napp, ... (n_app/e_app)
      sparse(N.node_app,N.node_app)]; %(n_app/n_app)
b = [ Us; ... (e_obj)
      -P*Is
      Us_app; ... (e_app)
      -(1j*omega2)*Is_napp]; %(n_app)  
toc
disp(' ')
%% solve
disp('Solving system ...')
tic
x = M\b;
toc
disp(' ')
%% extract solution
sol.I_obj = x(1:nSticks); % currents in the stick
sol.I_app = x(nSticks+nNodes+1:nSticks+nNodes+N.lumped); % currents in the lumped branches
sol.Q = (Aobj*sol.I_obj+A_eapp*sol.I_app+Is)/(1j*omega); % charges 
sol.U_obj = x(nSticks+1:nSticks+nNodes); % electric scalar potential in the nodes
sol.U_app = x(nSticks+nNodes+N.lumped+1:end);% electric scalar potential in the lumped nodes
%%
figure
subplot(2,1,1)
plot(real(sol.U_obj),'r o-')
hold on
plot(imag(sol.U_obj),'b o-')
title('Potential [V]')
subplot(2,1,2)
plot(real(sol.I_obj),'r o-')
hold on
plot(imag(sol.I_obj),'b o-')
title('Current [A]')
drawnow
%%
%% post processing
disp('Post-processing ...')
tic
%current density
Jre = zeros(3,nSticks);
Jim = zeros(3,nSticks);
for k=1:nSticks
    ut= (NN(:,G(2,k))-NN(:,G(1,k)))/lgt(k);
    Jre(:,k) = real(sol.I_obj(k))*ut./lgt(k);
    Jim(:,k) = imag(sol.I_obj(k))*ut/lgt(k);
end
toc
disp(' ')
%% plot
% J Re
cmax=max(max(real(sol.I_obj)));
cmin=min(min(real(sol.I_obj)));
figure
hold on
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color',[0.65 0.65 0.65])
quiver3_c(bar_stick(:,1).',bar_stick(:,2).',bar_stick(:,3).',(Jre(1,:)),(Jre(2,:)),(Jre(3,:)),real(sol.I_obj));
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
title('J vec \Re')
colorbar;
caxis([cmin cmax]);
view(3)
cmax=max(max(imag(sol.I_obj)));
cmin=min(min(imag(sol.I_obj)));
% J Im
figure
hold on
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color',[0.65 0.65 0.65])
quiver3_c(bar_stick(:,1).',bar_stick(:,2).',bar_stick(:,3).',(Jim(1,:)),(Jim(2,:)),(Jim(3,:)),imag(sol.I_obj));
axis equal
zlim auto
xlabel('x')
ylabel('y')
zlabel('z')
title('J vec imag')
colorbar;
caxis([cmin cmax]);
view(3)
% Potential Re
figure
cmin=min(min(real(sol.U_obj)));
cmax=max(max(real(sol.U_obj)));
scatter3(Matrix_P0(:,1),Matrix_P0(:,2),Matrix_P0(:,3),3,'CData',real(sol.U_obj),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title(['Phi_e \Re'])
axis([xmin xmax ymin ymax zmin zmax])
caxis([cmin cmax])
% Potential Im
figure
cmin=min(min(imag(sol.U_obj)));
cmax=max(max(imag(sol.U_obj)));
scatter3(Matrix_P0(:,1),Matrix_P0(:,2),Matrix_P0(:,3),3,'CData',imag(sol.U_obj),'MarkerFaceColor','flat')
axis equal
colorbar
colormap jet
xlabel('x')
ylabel('y')
zlabel('z')
view(3)
title(['Phi_e \Im'])
caxis([cmin cmax])
axis([xmin xmax ymin ymax zmin zmax])
drawnow
%% scattered and total fields in target points 
if E_target
if E_external==1
    sol.Eext=E_ext(TP(:,1),TP(:,2),TP(:,3));
else
    sol.Eext=zeros(size(TP,1),3);
end
[A] = fun_A_compute_vec(nSticks,TP,G.',NN.',sol.I_obj,0,freq);
[mgradphi] = fun_gradphi_compute_vec_ret(TP,NN.',sol.Q,freq);
sol.Esca=-1j*omega*A+mgradphi;
sol.Etot=sol.Esca+sol.Eext;
%
NNTP=[NN,TP.'];
xmin=min(NNTP(1,:))-1e-3;ymin=min(NNTP(2,:))-1e-3;zmin=min(NNTP(3,:))-1e-3;xmax=max(NNTP(1,:))+1e-3;ymax=max(NNTP(2,:))+1e-3;zmax=max(NNTP(3,:))+1e-3;
%
figure
hold on
plot3(TP(:,1),TP(:,2),TP(:,3),'.','markersize',15)
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
eeR = real(sol.Etot);
normER=sqrt(eeR(:,1).^2+eeR(:,2).^2+eeR(:,3).^2);
cmin=min(normER);
cmax=max(normER);
quiver3_c(TP(:,1),TP(:,2),TP(:,3),eeR(:,1),eeR(:,2),eeR(:,3),normER);
title('sol.Etot vec \Re  [V/m]')
axis equal
caxis([cmin cmax])
colorbar
axis([xmin xmax ymin ymax zmin zmax])
view(3)
figure
hold on
plot3(TP(:,1),TP(:,2),TP(:,3),'.','markersize',15)
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
eeI = imag(sol.Etot);
normEI=sqrt(eeI(:,1).^2+eeI(:,2).^2+eeI(:,3).^2);
cmin=min(normEI);
cmax=max(normEI);
quiver3_c(TP(:,1),TP(:,2),TP(:,3),eeI(:,1),eeI(:,2),eeI(:,3),(normEI));
title('sol.Etot vec \Im [V/m]')
axis equal
caxis([cmin cmax])
colorbar
view(3)
axis([xmin xmax ymin ymax zmin zmax])
end
%% show J figure
warning off
figure(2)
warning on
%% delete log files
delete log_dist.txt log_dist_ee.txt logL.txt logP.txt
%%

                           
                           
