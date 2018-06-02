addpath(genpath('C:\Users\thano\OneDrive\Documents\USC\DeepSolar'))
%% block matrix with ones(m)
clc; clear all; close all;
%m # of MESS, d # of days
mg = 4; days = 6;
MESS = 2;
if MESS > mg
    error('No of MESS > of Micro-grids')
end
% create theadjacency matrix Ad_g
%1st: create block diagonal matrix with ones sub-blocks
Ad_g = kron(eye(days-1),ones(mg));
%2nd: include zeros in the 1st m columns of the matrix
Ad_g = [zeros(mg*days-mg,mg) Ad_g];
%3ed: include zeros in the last m rows of the matrix
Ad_g = [Ad_g; zeros(mg,mg*days)];
% G is the graph
G = digraph(Ad_g);
% assign random weigths to the edges
% G.Edges.Weight = randi(4,numedges(G),1);
% add two nodes S* and T*
G = addnode(G,{'S*'});
G = addnode(G,{'T*'});
% G = addnode(G,2);
% add edges from S* to the 1st columns of nodes 
% and last column of nodes to T*
% 1st create the table with new edges
Ed_S = table([(mg*days+1)*ones(1,mg);1:mg ]',ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T = table([(mg*(days-1)+1):(mg*days) ;(mg*days+2)*ones(1,mg)]',ones(mg,1),'VariableNames',{'EndNodes','Weight'});
% 2nd add edges to graph
G = addedge(G,Ed_S);
G = addedge(G,Ed_T);
% connect S* and T* (add edge from m*d+2 to m*d+1 node
% G = addedge(G,mg*days+2,mg*days+1,MESS);
% assign random weigths to the edges
% G.Edges.Weight(1:mg*(days-1)) = randi(4,mg*(days-1),1);

%% add capacities and costs to each edge
% add Edge Label
G.Edges.Label = [1:numedges(G)]';
% init. capacities
G.Edges.Capacities = zeros(numedges(G),1);
% capacities if S* and T* are connected
% G.Edges.Capacities = [ ones(mg*(days-1)*mg,1) ; ones(numedges(G)- mg*mg*(days-1)-1,1); MESS];
% capacites with supply and demand in S* and T*
G.Edges.Capacities = [ ones(mg*(days-1)*mg,1) ; ones(numedges(G)- mg*mg*(days-1),1)];
% init. costs
G.Edges.Costs = zeros(numedges(G),1);
% edge cost are zero in edges connecting S* and T*
G.Edges.Costs = [randi(4,mg*(days-1)*mg,1) ; zeros(numedges(G)- mg*mg*(days-1),1)];
G.Edges.Weight = G.Edges.Costs;
% plot the graph
figure(1)
% set(gcf, 'Position', get(0, 'Screensize'));
h1 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h1,'layered','Direction','right','Sources',mg*days+1,'Sinks',mg*days+2);
%%
suc_sh_pa = zeros(1,MESS);
figure(2)
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(G,'S*','T*');
    suc_sh_pa(i_rm) = path_len;
    % set path cost to inf, to chose 2nd shortest path
    G.Edges.Weight(path1) = inf;
    highlight(h2,'Edges',path1,'EdgeColor','r')
    figure(i_rm+2)
    h2 = plot(G,'EdgeLabel',G.Edges.Weight);
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
    % print path
    disp(path1)
end
suc_sh_pa;
disp(['Total cost is  ',num2str(sum(suc_sh_pa))])
%%
suc_sh_pa = zeros(1,MESS);
figure(2)
h2 = plot(G,'EdgeLabel',G.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
% for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(G,'S*','T*');
%     suc_sh_pa(i_rm) = path_len;
    Rm_nodes = P_nodes(2:(end-1));
%     G = rmnode(G,Rm_nodes);
    highlight(h2,'Edges',path1,'EdgeColor','r')
%     figure(i_rm+2)
    h2 = plot(G,'EdgeLabel',G.Edges.Weight);
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
% end
suc_sh_pa;
disp(['Total cost is  ',num2str(sum(suc_sh_pa))])
%%
close all;
%% ~~~~~~~~ Test codes ~~~~~~~~~~~~~~~~~~

%%
% Test 2
% kron(eye(days-1),ones(mg));
% % kron(eye(2),ones(2))
% eye(mg) is the diagonal with the cost difference with & w/o MESS
% ones is the matrix denoting the tranfer cost from microgrid i to j
A0 = [eye(mg) zeros(mg);zeros(mg) ones(mg)];
AD0 = kron(eye(days-1),A0);
% spy(AD0)
AD1 = blkdiag(AD0,eye(mg));
AD2 = [zeros((days-1)*2*mg+mg,mg) AD1; zeros(mg) zeros(mg,(days-1)*2*mg+mg)];
% spy(AD2)
G0 = digraph(AD2);
G0 = addnode(G0,{'S*'});
G0 = addnode(G0,{'T*'});
Ed_S0 = table([(2*mg*days+1)*ones(mg,1) (1:mg)'],mg*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T0 = table([ ((2*mg*days-mg+1):(2*mg*days))' (2*mg*days+2)*ones(mg,1)],mg*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
G0 = addedge(G0,Ed_S0);
G0 = addedge(G0,Ed_T0);
figure(232)
h3 = plot(G0);
layout(h3,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% highlight(h3,'Edges',,'EdgeColor','r')
%% Test 1
B1 = [eye(2); zeros(6,2)];
B2 = [zeros(2); eye(2); zeros(4,2)];
Ad_g0 = full(adjacency(G));
Ad_g1 = [Ad_g0(:,1:2) B1 Ad_g0(:,3:4) B2 Ad_g0(:,5:8) ]
G0 = digraph(Ad_g1);
figure(124)
h12 = plot(G0)
layout(h12,'layered','Direction','right')
%%
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