addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
%%
addpath(genpath('C:\Users\thano\OneDrive\Documents\USC'))
cd C:\Users\thano\OneDrive\Documents\USC\DeepSolar\sustech\MESS-optimization-master\MESS-optimization
%%
% graph with node split
clc; clear all; close all;
%m , d # of days
mg = 3; days = 5;
% # of MESS
MESS = 2;
if MESS > mg
       error('No of MESS > of Micro-grids')
end
% eye(mg) is the diagonal with the cost difference with & w/o MESS
% ones is the matrix denoting the tranfer cost from microgrid i to j
A0 = [eye(mg) zeros(mg);zeros(mg) ones(mg)];
AD0 = kron(eye(days-1),A0);
% spy(AD0)
AD1 = blkdiag(AD0,eye(mg));
AD2 = [zeros((days-1)*2*mg+mg,mg) AD1; zeros(mg) zeros(mg,(days-1)*2*mg+mg)];
% spy(AD2)
G0 = digraph(AD2);
% initial number of edges
init_edges = numedges(G0);
G0 = addnode(G0,{'S*'});
G0 = addnode(G0,{'T*'});
Ed_S0 = table([(2*mg*days+1)*ones(mg,1) (1:mg)'],MESS*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T0 = table([ ((2*mg*days-mg+1):(2*mg*days))' (2*mg*days+2)*ones(mg,1)],MESS*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
G0 = addedge(G0,Ed_S0);
G0 = addedge(G0,Ed_T0);
%% add edges costs and capacities 
% edge capacities are 1 except those of S* and T* with capacity equal to
% the # of MESS
G0.Edges.Capacities = [ones(init_edges,1); MESS*ones(length((init_edges+1):numedges(G0)),1) ];
% edges have random cost; edges from S* and T* have zero costs
G0.Edges.Costs = [ randi(4,init_edges,1) ; zeros(length((init_edges+1):numedges(G0)),1)];
% G0.Edges.Costs = [ 4*rand(init_edges,1) ; zeros(length((init_edges+1):numedges(G0)),1)];
G0.Edges.Weight = G0.Edges.Costs;
% assign label nodes
G0.Edges.Labels = (1:numedges(G0))';
figure(1)
h1 =plot(G0,'EdgeLabel',G0.Edges.Costs);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% highlight(h3,'Edges',,'EdgeColor','r')
%%
cost_v = G0.Edges.Costs';
tic;
cvx_begin
cvx_solver gurobi
variable fl(numedges(G0),1)  %integer
% maximize sum(cost_v*fl)
minimize cost_v*fl
subject to 
    incidence(G0)*fl == [zeros(numnodes(G0)-2,1); -MESS; MESS];
    0 <= fl <= G0.Edges.Capacities;
cvx_end
toc;
[fl cost_v'];
disp(['Total cost: ',num2str(cost_v*fl)])
%%
highlight(h1,'Edges',find(fl>0),'EdgeColor','r','LineWidth',1.5)
%% find path of nodes 
% find start and end nodes of edges with flow
[st_ed,end_ed] = findedge(G0);
% create matrix of nodes
ed_mat1 = [st_ed end_ed];
% remove the last 2*mg rows which include nodes S* and T*
ed_mat = ed_mat1(1:(end-2*mg),:);
fl1 = fl(1:(end-2*mg),:);
% index matrix of the nodes of the graph with flow
ind_ed = ed_mat((fl1>0),:);
% reshape the matrix to chose every other row
res_ind_ed = [reshape(ind_ed',2*MESS,[])]';
% chooose even rows of the reshaped array until no of rows size(res_ind_ed,1)
res_ind_ed = res_ind_ed(2:2:size(res_ind_ed,1),:);
% reshape the even rows of the reshaped matrix
ind_ed2 = [reshape(res_ind_ed',2,[])]';
%% mod ind_ed with mg to find which micro-grid it belonds
ind_mod_ed = mod(ind_ed2,mg);
% if any element is 0 then it is the mg-th micro-grid
ind_mod_ed(ind_mod_ed == 0) = mg;
% final matrix of micro-grid scheduling
mic_mat = ind_mod_ed;
%% print result from LP min cost 
clc;
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Day 1')
for i_mic = 1:(days-1)*MESS
    disp(['Transfer from ',num2str(mic_mat(i_mic,1)),' to ',num2str(mic_mat(i_mic,2))]);
    if mod(i_mic,MESS) == 0 && i_mic ~= (days-1)*MESS
        disp(['Day ',num2str((i_mic/MESS)+1)]);
    end
end
disp(['Total cost with LP is: ',num2str(G0.Edges.Costs'*fl)])

%%
Gs = G0;
suc_sh_pa = zeros(1,MESS);
figure(2)
h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*');
    suc_sh_pa(i_rm) = path_len;
    Rm_nodes = P_nodes(2:(end-1));
%     G = rmnode(G,Rm_nodes);
    Gs.Edges.Weight(path1) = inf;
    highlight(h2,'Edges',path1,'EdgeColor','r','LineWidth',1.5)
%     figure(i_rm+2)
%     h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
end
suc_sh_pa;
disp(['Total cost with shortesth paths is  ',num2str(sum(suc_sh_pa))])
%% find the paths from LP min cost flow solution 
for i_path = 1:MESS
    mic_mat(i_path:MESS:size(mic_mat,1),:)
end
%%
clc
s_path = Gs.Edges.Labels(fl>0);
[psOut,ptOut] = findedge(Gs,s_path);
[psOut ptOut];
for i_n_path = 1:length(psOut)
    next_node_index = find(ptOut(i_n_path) == psOut(:,:))
    [psOut(next_node_index) ptOut(next_node_index)] 
    i_n_path = next_node_index;
end

%%
 ind_ed_s = ed_mat((fl1>0),:);
%  number of edges between S* and T* w/o weights
 dis_S_T = distances(Gs,numnodes(Gs)-1,numnodes(Gs),'Method','unweighted');
 % distance between day 1 and last day:
 dis_d1_dl = dis_S_T - 2;
for i_MESS = 1:MESS
    disp(['+++++++++ Path #',num2str(i_MESS),'++++++++++++++'])
% the next path is the st_ind row:
% st_ind = find(ind_ed_s(:,1) == ind_ed_s(i_MESS,2))
   
    st_ind = find(ind_ed_s(:,1) == ind_ed_s(1,2));
    np = [0 0];
%     np = ind_ed_s(st_ind,:)
    i_path = 0;
    for i__path = 1:(dis_d1_dl-1)
        np = ind_ed_s(st_ind,:);
        % next node is: 
        nex_nd = np(2);
        % iterate
        st_ind = find(ind_ed_s(:,1) == nex_nd);
        i_path = i_path + 1;
        % path_mat(:,i_path) = np
        path_mat(i_path,:) = np;
    end
    fin_path_mat = [ind_ed_s(1,:) ; path_mat]
    % node ID's of the 1st path
    node_ids = unique(fin_path_mat);
    % end
    % find edges with start & end nodes in path mat"
    edges_id = findedge(Gs,fin_path_mat(:,1),fin_path_mat(:,2));
    path_cost = sum(Gs.Edges.Costs(edges_id));
    % remove the node ID's from the original matrix and repeat
    % compare the difference in the elements of each column individualy
    dif_col1 = setdiff(ind_ed_s(:,1),fin_path_mat(:,1));
    dif_col2 = setdiff(ind_ed_s(:,2),fin_path_mat(:,2));
    dif_col = [dif_col1 dif_col2];
    %dif col is the new ind_ed_s matrix
   ind_ed_s = dif_col;
   Aaaasdads = 523424;
end
%% plot each path with green in new figure
figure(3)
h3 = plot(Gs);%,'EdgeLabel',Gs.Edges.Weight);
layout(h3,'layered','Direction','right','Sources','S*','Sinks','T*')
highlight(h3,'Edges',edges_id,'EdgeColor','g','LineWidth',1.5)
title('Plot ind. paths')

%% extract sub-graph from G0 graph with the shortest paths 
% ind_ed_s =  ed_mat1(fl>0,:);
% H0 = subgraph(G0,unique(ind_ed_s));
H0 = subgraph(G0,unique(path_mat));
H0 = addnode(H0,{'S*'});
% H0 = addnode(H0,{'T*'});
figure(443)
h55 =plot(H0,'EdgeLabel',H0.Edges.Costs);
layout(h55,'layered','Direction','right','Sources', 'S*','Sinks','T*')
[diff_paths test_car]= conncomp(H0,'Type','weak')


