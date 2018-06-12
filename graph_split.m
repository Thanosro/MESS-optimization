addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
rmpath('C:\Users\Thanos\Documents\DeepSolar\Optimal_flow\cvx\lib\narginchk_')
%%
addpath(genpath('C:\Users\thano\OneDrive\Documents\USC'))
cd C:\Users\thano\OneDrive\Documents\USC\DeepSolar\sustech\MESS-optimization-master\MESS-optimization
rmpath('C:\Users\thano\OneDrive\Documents\USC\DeepSolar\OPF\cvx\lib\narginchk_')
%%
% graph with node split
clc; clear all; close all;
%m , d # of days
mg = 4; days =7;
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
title('Min cost flow LP')
% highlight(h3,'Edges',,'EdgeColor','r')
%%
cost_v = G0.Edges.Costs';
tic;
cvx_begin quiet
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
disp(['Total cost with LP: ',num2str(cost_v*fl)])
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
% clc;
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
%% DEBUGGING 
% compare edges costs
isequal(Gs.Edges.Weight,G0.Edges.Costs)
%% SHORTEST PATHS
Gs = G0;
suc_sh_pa = zeros(1,MESS);
figure(2)
h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
title('Shortest Paths')
disp('******* Shortest Paths Costs ***********')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*');
    suc_sh_pa(i_rm) = path_len;
    Rm_nodes = P_nodes(2:(end-1));
%     G = rmnode(G,Rm_nodes);
    Gs.Edges.Weight(path1) = inf;
    highlight(h2,'Edges',path1,'EdgeColor','r','LineWidth',1.5)
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
    disp(['Cost of ',num2str(i_rm),' path is: ',num2str(suc_sh_pa(i_rm))])
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % each column of Shortest_Path_mat is the edges index for each path 
    Shortest_path_mat(:,i_rm) = path1;
end
suc_sh_pa;
disp(['Total cost with shortest paths is  ',num2str(sum(suc_sh_pa))])
disp(newline)
%% DEBUGGING check if the paths are the same (LP & shortest paths)
disp('**************** Debug Section ************')
if isequal(fl,isinf(Gs.Edges.Weight))== 1
    disp('Paths are the same')
else
    disp('Paths are not the same')
end
if sum(suc_sh_pa) == (G0.Edges.Costs'*fl) 
    disp('Costs are equal')
else
    disp('Costs are not equal')
end
disp('******************************')
disp(newline)
%% find the paths from LP min cost flow solution 
ind_ed_s = ed_mat((fl1>0),:);
%  number of edges between S* and T* w/o weights
dis_S_T = distances(G0,numnodes(G0)-1,numnodes(G0),'Method','unweighted');
 % distance between day 1 and last day:
 dis_d1_dl = dis_S_T - 2;
%  edges_id_mat = zeros(1,MESS);
 disp('******* LP min cost individual Paths ***********')
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
    fin_path_mat = [ind_ed_s(1,:) ; path_mat];
    % node ID's of the 1st path
    node_ids = unique(fin_path_mat);
    % end
    % find edges with start & end nodes in path mat"
    edges_id = findedge(G0,fin_path_mat(:,1),fin_path_mat(:,2));
    path_cost(i_MESS) = sum(G0.Edges.Costs(edges_id)); 
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % each column of LP_path_mat is the edges index for each path 
    LP_path_mat(:,i_MESS) = edges_id';
    disp(['Cost of ',num2str(i_MESS),' path is: ',num2str(path_cost(i_MESS))])
    % remove the node ID's from the original matrix and repeat
    % compare the difference in the elements of each column individualy
    dif_col1 = setdiff(ind_ed_s(:,1),fin_path_mat(:,1));
    dif_col2 = setdiff(ind_ed_s(:,2),fin_path_mat(:,2));
    dif_col = [dif_col1 dif_col2];
    %dif col is the new ind_ed_s matrix
   ind_ed_s = dif_col;
end
disp(['Total cost is with LP: ',num2str(sum(path_cost))])
disp(newline)
%% DEBUG compare pahts from 2 solutions
% remove 1st and last row of Shortest path mat cuz it's the edges S* T*
% if condition holds to run the code more than once
if size(Shortest_path_mat,1) ~= size(LP_path_mat,1)
    Shortest_path_mat(1,:) = [];
    Shortest_path_mat(end,:) = [];
end
% edge index of each approach
% edge indices of each path from LP is LP_path_mat
% edge indices of each path from Shortest paths is Shortest_path_mat
% [LP_path_mat Shortest_path_mat]
% cost of each edge of each path 
disp('Cost of each edge of LP is:')
disp(num2str(G0.Edges.Costs(LP_path_mat)))
% G0.Edges.Costs(LP_path_mat)
disp('Cost of each edge with Shortest Paths is:')
disp(num2str(Gs.Edges.Costs(Shortest_path_mat)))
% Gs.Edges.Costs(Shortest_path_mat)
% cost of each path with each approach
disp(['Cost of each path of LP is: ',newline,num2str(sum(G0.Edges.Costs(LP_path_mat)))])
% sum(G0.Edges.Costs(LP_path_mat))
disp(['Cost of each path with Shortest paths is:',newline,num2str(sum(Gs.Edges.Costs(Shortest_path_mat)))])
% sum(Gs.Edges.Costs(Shortest_path_mat))
disp('The increasing cost of each path with LP is:')
disp(num2str(cumsum(G0.Edges.Costs(LP_path_mat))))
% cumsum(G0.Edges.Costs(LP_path_mat))
disp('The increasing cost of each path with Shortest Path is')
disp(num2str(cumsum(Gs.Edges.Costs(Shortest_path_mat))))
% cumsum(Gs.Edges.Costs(Shortest_path_mat))
%% DEBUGGING plot the increasing cost 
for i_plot = 1:MESS
    figure(46+i_plot)
    plot([cumsum(Gs.Edges.Costs(Shortest_path_mat(:,i_plot))) cumsum(G0.Edges.Costs(LP_path_mat(:,i_plot))) ])
    title(['Increasing cost of path # ',num2str(i_plot)])
    xlabel('# of edges')
    ylabel('Cost')
    legend('SP','LP')
end
figure(100)
plot([sum([cumsum(G0.Edges.Costs(LP_path_mat))],2) sum(cumsum(Gs.Edges.Costs(Shortest_path_mat)),2)])
title('Comparison increasing path cost')
xlabel('# of edges')
ylabel('Cost')
legend('LP','SP')
%% DEGUBBING FIND FROM EACH NODE NEXT AVAILABLE NODES & EDGES' COSTS
% LP %
ind_ed2(:,1) % shows nodes in each path chosen by LP
for i_nx_nd = 1:2:size(ind_ed2,1)
    disp('Index of available edges is:')
    outedges(G0,ind_ed2(i_nx_nd,1))
    disp('Cost of available edges is:')
    G0.Edges.Costs(outedges(G0,ind_ed2(i_nx_nd,1)))
    % print which edge the algorithm chose 
%     disp(['LP chose edges # ',num2str(),' with cost',num2str()])
end
% G0.Edges.Costs(outedges(G0,ind_ed(i_nx_nd,1)))
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
%% DEBUG check for different edges
% test the different edges in the 2 approaches
% edges in shortest paths: isinf(Gs.Edges.Weight)
% edges in in LP min-cost flow: fl
dif_egd = [isinf(Gs.Edges.Weight) fl]
% compare between the 2 columns and return the index of the differences
(dif_egd(:,1) == dif_egd(:,2)) == 0
(dif_egd(:,1) == 1) == (dif_egd(:,2) == 1)

[find(dif_egd(:,1) == 1) find(dif_egd(:,2) == 1)]


%% test code
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
%% test code edge label 
Ed_lbls = reshape(G0.Edges.Labels,mg,[]);
% add 1 column at start
Ed_lbls = [zeros(mg,1) Ed_lbls(:,:)];
Ed_lbls(:,1) = Ed_lbls(:,end);
Ed_lbls(:,end) = [];
Ed_lbls =  flipud(Ed_lbls);
%% matrix with zeros same as Ed_lbls to show the paths
Ind_G0 = zeros(size(Ed_lbls));
% Ind_G0(fl>0) = 1;
Ind_Gs = zeros(size(Ed_lbls));

%%
highlight(h1,'Edges',68,'EdgeColor','r','LineWidth',1.5)
