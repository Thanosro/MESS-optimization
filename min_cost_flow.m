% syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 real 
% x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]';
% x = randi(4,numedges(G),1);
% incidence(G)*x
% %%
% cost_v = randi(5,1,numedges(G));
%%
cost_vec = G.Edges.Costs';
cvx_begin
cvx_solver gurobi
variable fl(numedges(G),1)

% maximize sum(cost_v*fl)
minimize cost_vec*fl
subject to 
    incidence(G)*fl == [ zeros((numnodes(G)-2),1); -MESS; MESS];
    0 <= fl <= G.Edges.Capacities;
cvx_end
[fl cost_vec'];
 %% show edges on graph plot
highlight(h1,'Edges',find(fl>0),'EdgeColor','r','LineWidth',1.5)
%% find path of nodes 
% find start and end nodes of edges with flow
[st_ed,end_ed] = findedge(G);
% create matrix of nodes
ed_mat1 = [st_ed end_ed];
% remove the last 2*mg rows which include nodes S* and T*
ed_mat = ed_mat1(1:(end-2*mg),:);
fl1 = fl(1:(end-2*mg),:);
% index matrix of the nodes of the graph with flow
ind_ed = ed_mat((fl1>0),:);
% mod ind_ed with mg to find which micro-grid it belonds
ind_mod_ed = mod(ind_ed,mg);
% if any element is 0 then it is the mg-th micro-grid
ind_mod_ed(ind_mod_ed == 0) = mg;
% final matrix of micro-grid scheduling
mic_mat = ind_mod_ed;
%% print result 
clc;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Day 1')
for i_mic = 1:(days-1)*MESS
    disp(['Transfer from ',num2str(mic_mat(i_mic,1)),' to ',num2str(mic_mat(i_mic,2))]);
    if mod(i_mic,MESS) == 0 && i_mic ~= (days-1)*MESS
        disp(['Day ',num2str((i_mic/MESS)+1)]);
    end
end
disp(['Total cost is ',num2str(G.Edges.Costs'*fl)])
%%
% reproduce result with:
mg = 7; days = 4;
MESS = 4;
G.Edges.Costs = [1;3;1;1;3;1;4;4;3;1;3;3;4;3;4;2;2;4;1;1;...
    1;2;4;4;1;2;3;2;3;3;2;2;1;4;1;1;2;1;2;2;4;4;1;3;2;2;3;4;2;...
   4;2;3;3;3;3;3;1;1;4;1;1;3;4;3;1;2;2;4;1;4;3;2;1;2;2;1;3;1;2;...
    3;2;2;3;2;4;4;3;2;3;1;4;4;4;2;3;1;2;2;1;1;2;1;3;2;3;3;3;1;1;2;3;3;...
    2;4;3;4;3;2;1;3;4;2;1;2;1;2;2;3;2;4;3;4;3;4;1;3;2;3;3;1;2;1;3;4;2;4;3;...
    0;0;0;0;0;0;0;0;0;0;0;0;0;0];
%% draft
flow_nodes = G.Edges.EndNodes([fl>0]);
path_vec = findnode(G,flow_nodes);
path_list = mod(path_vec,mg)
reshape(path_vec,MESS,[])
%%
ind_mat = 1:mg*days;
ind_mat = reshape(ind_mat,mg,days);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
for i_mess = 1:length(path_vec)-2
    [mg_no, day_no] = find(ind_mat == path_vec(i_mess)); 
    disp(['Day ',num2str(day_no),' at Micro-grid ',num2str(mg_no)])
end

