function [NN,G,nNodes,nSticks,Matrix_G] = fun_extract_from_comsol_1D(C)
disp('reading sources'' geo data from comsol txt...')
tic
for ii = 1:size(C,1) 
    CC = char(C(ii));
   if size(CC) == size('# Mesh point coordinates')%count number of characters
       if CC ==  '# Mesh point coordinates'
           ind_start_point = ii+1;
       end
   end
   if size(CC) == size('2 # number of nodes per element')
       if CC ==  '2 # number of nodes per element'
           ind_start_stick = ii+3;
       end
   end
end
%%
CC=char(C(19));
ind=find(CC=='#');
nNodes=str2num(CC(1:ind-1));
%%
CC=char(C(ind_start_stick-2));
ind=find(CC=='#');
nSticks=str2num(CC(1:ind-1));
%%
NN = str2num(char(C(ind_start_point:ind_start_point+nNodes)))';
G = str2num(char(C(ind_start_stick:ind_start_stick+nSticks)))'; 
G = G+1; %comsol starting node = 0 ---> becomes 1
%%
%compute full G matrix
val_neg_g=G(1,:);
val_pos_g=G(2,:);
Matrix_G=sparse([1:nSticks],val_neg_g,-ones(nSticks,1),nSticks,nNodes);
Matrix_G=Matrix_G+sparse([1:nSticks],val_pos_g,ones(nSticks,1),nSticks,nNodes);
%%
check_node = max(max(G)) == nNodes; %if check = 1 OK
if check_node == 0
    disp('ERROR: the number of nodes does not match')
    return
end
toc
disp('done')
disp('---------------------------------')
%%
disp('mesh has:')
disp([num2str(nSticks) ' edge elements'])
disp('---------------------------------')
%% plot
disp('computing baricentres and plotting...')
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
% tic
% bar_stick=zeros(3,nSticks);
% for k=1:nSticks
%     bar_stick(:,k)=(NN(:,G(1,k))+NN(:,G(2,k)))/2;
% end
% toc

end

