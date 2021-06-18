clear
close all
clc
restoredefaultpath
%%
Lumped=[];Lumped.node_start=[];Lumped.node_end=[];Lumped.R=[];Lumped.L=[];Lumped.Cinv=[];Lumped.value=[];
Current_Source=[];Current_Source.node=[];Current_Source.value=[];
%% BEGIN USER SETTINGS
% generate below variable NN (3,nNodes) which contains the coordinate of the stick mesh
long=1; % length of the two_wire_line [m]
dist=0.2; % distance between the two lines [m]
nsp=100; % number of points
NN1=zeros(3,nsp);
NN1(1,:)=linspace(0,long,nsp);
NN1(3,:)=dist*0.5;
NN2=zeros(3,nsp);
NN2(1,:)=linspace(0,long,nsp);
NN2(3,:)=-dist*0.5;
NN=[NN1,NN2];
% generate below variable GG (2,nSticks) which is the edges-nodes incidence matrix
G=[1:nsp-1;2:nsp];
G=[G,G+nsp];
% generate Lumped element (R, L, C, or voltage source elements)
% Lumped(1).node_start=-1; % node1 of the lumped branch (if positive it is a mesh node, if negative is an appended node)
% Lumped(1).node_end=1;  % node2 of the lumped branch (if positive it is a mesh node, if negative is an appended node)
% Lumped(1).R=0; % lumped resitance of the lumped branch
% Lumped(1).L=0; % lumped inductance of the lumped branch
% Lumped(1).Cinv=0; % lumped "inverse of  the capacitance" of the lumped branch
% Lumped(1).value=1; % voltage source value 
% Lumped(2).node_start=-1;
% Lumped(2).node_end=nsp+1;
% Lumped(2).R=0;
% Lumped(2).L=0;
% Lumped(2).Cinv=0;
% generate injected currents
Current_Source(1).node=1; % node where the current is injected 
Current_Source(1).value=1; % value of the injected current
Current_Source(2).node=101;
Current_Source(2).value=-1;
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
%%

