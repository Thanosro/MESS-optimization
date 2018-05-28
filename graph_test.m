% block matrix with ones(m)
clc; clear all; close all;
%m # of MESS, d # of days
m = 4; d = 5;
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
Ed_S = table([(m*d+1)*ones(1,m);1:m ]',ones(m,1),'VariableNames',{'EndNodes','Weight'});
Ed_T = table([(m*(d-1)+1):(m*d) ;(m*d+2)*ones(1,m)]',ones(m,1),'VariableNames',{'EndNodes','Weight'});
% 2nd add edges to graph
G = addedge(G,Ed_S);
G = addedge(G,Ed_T);
% connect S* and T* (add edge from m*d+2 to m*d+1 node
G = addedge(G,m*d+2,m*d+1,1);
% plot the graph
figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
h = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h,'layered','Direction','right','Sources',m*d+1,'Sinks',[m*d+2])
%%
spy(adjacency(G))
%%
spy(incidence(G))