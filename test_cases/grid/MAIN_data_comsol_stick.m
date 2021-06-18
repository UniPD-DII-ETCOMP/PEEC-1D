%% ****************************** data_comsol *****************************
clear
close all
clc
%%
Lumped=[];Lumped.node_start=[];Lumped.node_end=[];Lumped.R=[];Lumped.L=[];Lumped.Cinv=[];Lumped.value=[];
Current_Source=[];Current_Source.node=[];Current_Source.value=[];
%% BEGIN USER SETTINGS
% generate variable NN (3,nNodes) which contains the coordinate of the stick mesh
C = textread('data_com.mphtxt', '%s','delimiter', '\n'); % load data from comsol 
[NN,G,~,~,~] = fun_extract_from_comsol_1D(C);%main_extract;
% generate variable GG (2,nSticks) which is the edges-nodes incidence matrix
% already done above
% generate Lumped element (R, L, C, or voltage source elements)
Lumped(1).node_start=1; % corner 1
Lumped(1).node_end=1001; % corner 3
Lumped(1).R=0;
Lumped(1).L=0;
Lumped(1).Cinv=0;
Lumped(1).value=1;
% generate injected currents
% Current_Source(1).node=1; % node where the current is injected 
% Current_Source(1).value=1; % value of the injected current
% Current_Source(2).node=1001;
% Current_Source(2).value=-1;
% END USER SETTINGS
%%
nNodes=size(NN,2);
nSticks=size(G,2);
%%
figure(1)
plot3(0,0,0,'r')
hold on
plot3(0,0,0,'g')
plot3(0,0,0,'b')
axis equal
grid on
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
plot3([NN(1,G(1,:));NN(1,G(2,:))],[NN(2,G(1,:));NN(2,G(2,:))],[NN(3,G(1,:));NN(3,G(2,:))],'.-','color','k')
%% 
save data.mat G NN Lumped Current_Source







