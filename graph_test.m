% block matrix with ones(m)
clc; clear all; %close all;
%m # of MESS, d # of days
m = 3; d = 2;
% create theadjacency matrix Ad_g
%1st: create block diagonal matrix with ones sub-blocks
Ad_g = kron(eye(d-1),ones(m));
%2nd: include zeros in the 1st m columns of the matrix
Ad_g = [zeros(m*d-m,m) Ad_g];
%3ed: include zeros in the last m rows of the matrix
Ad_g = [Ad_g; zeros(m,m*d)];
% G is the graph
G = digraph(Ad_g);
% add two nodes S* and T*
% G = addnode(G,{'S*'})
% G = addnode(G,{'T*'})
G = addnode(G,2);
% add edges from S* to the 1st columns of nodes 
% and last column of nodes to T*
% 1st create the table with new edges
New_ed = table([4 5 6 ;8 8 8]',[1 1 1]','VariableNames',{'EndNodes','Weight'});
% 2nd add edges to graph
G = addedge(G,New_ed);
% plot the graph
figure(2)
% set(gcf, 'Position', get(0, 'Screensize'));
plot(G)
%%
% add edge node indices
% m*d+2 is T* (m*d+2)*ones(1,m)
% m*d+1 is s* (m*d+1)*ones(1,m)
% 1st column is 1:m
% last column is (m*(d-1)+1):m*d
(m*d+2)*ones(1,m)

%%
spy(adjacency(G))