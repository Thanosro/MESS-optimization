% block matrix with ones(m)
clc; clear all; close all;
%m # of MESS, d # of days
mg = 5; days = 6;
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
G = addnode(G,{'S*'});
G = addnode(G,{'T*'});
% G = addnode(G,2);
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
h1 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h1,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2);
%%
MESS = 3;
if MESS > mg
    error('No of MESS > of Micro-grids')
end
figure(2)
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(G,'S*','T*');
    Rm_nodes = P_nodes(2:(end-1));
    G = rmnode(G,Rm_nodes);
    highlight(h2,'Edges',path1,'EdgeColor','r')
    figure(i_rm+2)
    h2 = plot(G,'EdgeLabel',G.Edges.Weight);
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
end
%% ~~~~~~~~ Test codes ~~~~~~~~~~~~~~~~~~
% matlab shortesh path 
% for i_rm = 1:3
[P_nodes,path_len,path1] = shortestpath(G,(mg*days+1),(mg*days+2));
% end
highlight(h1,'Edges',path1,'EdgeColor','r')
%%
figure(2)
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2)
% close all;
Rm_nodes = P_nodes(2:(end-1));
G = rmnode(G,Rm_nodes);
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
%%
spy(adjacency(G))
%%
spy(incidence(G))