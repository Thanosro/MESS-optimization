% block matrix with ones(m)
clc; clear all; close all;
%m # of MESS, d # of days
mg = 4; days = 5;
% create theadjacency matrix Ad_g
%1st: create block diagonal matrix with ones sub-blocks
Ad_g = kron(eye(days-1),ones(mg));
%2nd: include zeros in the 1st m columns of the matrix
Ad_g = [zeros(mg*days-mg,mg) Ad_g];
%3ed: include zeros in the last m rows of the matrix
Ad_g = [Ad_g; zeros(mg,mg*days)];
% G is the graph
G = digraph(Ad_g);
% add two nodes S* and T*
% G = addnode(G,{'S*'})
% G = addnode(G,{'T*'})
G = addnode(G,2);
% add edges from S* to the 1st columns of nodes 
% and last column of nodes to T*
% 1st create the table with new edges
Ed_S = table([(mg*days+1)*ones(1,mg);1:mg ]',mg*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T = table([(mg*(days-1)+1):(mg*days) ;(mg*days+2)*ones(1,mg)]',mg*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
% 2nd add edges to graph
G = addedge(G,Ed_S);
G = addedge(G,Ed_T);
% connect S* and T* (add edge from m*d+2 to m*d+1 node
G = addedge(G,mg*days+2,mg*days+1,mg);
% plot the graph
figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
h = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2)
%% matlab shortesh path
% for i_rm = 1:3
% [path1,path_len] = shortestpath(G,(mg*days+1),(mg*days+2));
P_nodes= shortestpath(G,(mg*days+1),(mg*days+2));
% end
% highlight(h,path1,'EdgeColor','r')
%%
figure(2)
h = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2)
close all;
% G = rmnode(G,P_nodes(2:end-1));
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2)
%%
spy(adjacency(G))
%%
spy(incidence(G))