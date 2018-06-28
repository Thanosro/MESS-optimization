% addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
% cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
% rmpath('C:\Users\Thanos\Documents\DeepSolar\Optimal_flow\cvx\lib\narginchk_')
% %% Laptop
% addpath(genpath('C:\Users\thano\OneDrive\Documents\USC'))
% cd C:\Users\thano\OneDrive\Documents\USC\DeepSolar\sustech\MESS-optimization-master\MESS-optimization
% rmpath('C:\Users\thano\OneDrive\Documents\USC\DeepSolar\OPF\cvx\lib\narginchk_')
%%
% graph with node split
% clc; clear all;
% clc;
clearvars -except MESS LP_cost ESS_cost_mat mg days
close all;
% load 100_mg_10d_Kwh.mat
%m , d # of days
% mg = 100; days =10;
% # of MESS
% MESS = 4;
if MESS > mg
       error('No of MESS > of Micro-grids')
end
% eye(mg) is the diagonal with the cost difference with & w/o MESS
% ones is the matrix denoting the tranfer cost from microgrid i to j
% create the relocation cost matrix
base_reloc_cost = 2.5;
reloc_mat =tril(ones(mg));
for dist = 1:mg
if dist > mg
    error('dist > mg')
end
reloc_mat(dist:mg+1:end) = dist;
end
reloc_mat = tril(reloc_mat);
reloc_mat =  base_reloc_cost*(reloc_mat+triu(reloc_mat',1));
A0 = [eye(mg) zeros(mg);zeros(mg) reloc_mat];
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
%%
% load 100_mg_10_days_MOND.mat
load 100_mg_10d_Kwh.mat
%%
% d_cost_red = -d_cost_red;
%% add edges costs and capacities 
% edge capacities are 1 except those of S* and T* with capacity equal to
% the # of MESS
% assign label nodes
G0.Edges.Labels = (1:numedges(G0))';
G0.Edges.Capacities = [ones(init_edges,1); MESS*ones(length((init_edges+1):numedges(G0)),1) ];
% edges have random cost; edges from S* and T* have zero costs
% G0.Edges.Costs = [ randi(4,init_edges,1) ; zeros(length((init_edges+1):numedges(G0)),1)];
% G0.Edges.Costs = [ 4*rand(init_edges,1) ; zeros(length((init_edges+1):numedges(G0)),1)];
G0.Edges.Costs = G0.Edges.Weight;
G0.Edges.Costs((numedges(G0)-2*mg+1):numedges(G0)) = 0;
% ----------------------------
hglht_edges = reshape(G0.Edges.Labels,mg,[]);
hglht_edges(:,(end-1:end)) =[];
% each column is the index daily cost of the matrix
dly_op_cost_ind = hglht_edges(:,(1:(mg+1):size(hglht_edges,2)));
dly_op_cost_ind = reshape(dly_op_cost_ind,[],1);
G0.Edges.Costs(dly_op_cost_ind) = d_cost_red;
% ----------------------------
G0.Edges.Costs(G0.Edges.Weight == base_reloc_cost) = 0;
%-----------------------------
G0.Edges.Weight = G0.Edges.Costs;
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
disp(['Total cost reduction with MESS: ',num2str(cost_v*fl)])
%%
% take the min of the average cost reduction of each micro-grid for 10 days
% avg_cost_red_mg = mean(reshape(d_cost_red,mg,days),2);
avg_cost_red_mg = sum(reshape(d_cost_red,mg,days),2);
[ESS_cost, MG_ESS_ind] = mink(avg_cost_red_mg,MESS);
disp(['Total cost reduction with ESS: ',num2str(sum(ESS_cost))])
% disp(['Total cost   with LP: ',num2str(cost_v*fl+Base_cost)])
%%
% highlight(h1,'Edges',find( fl>0),'EdgeColor','r','LineWidth',1.5)

% %% find path of nodes 
% % find start and end nodes of edges with flow
% [st_ed,end_ed] = findedge(G0);
% % create matrix of nodes
% ed_mat1 = [st_ed end_ed];
% % remove the last 2*mg rows which include nodes S* and T*
% ed_mat = ed_mat1(1:(end-2*mg),:);
% fl1 = fl(1:(end-2*mg),:);
% % index matrix of the nodes of the graph with flow
% ind_ed = ed_mat((fl1>0),:);
% % reshape the matrix to chose every other row
% res_ind_ed = [reshape(ind_ed',2*MESS,[])]';
% % chooose even rows of the reshaped array until no of rows size(res_ind_ed,1)
% res_ind_ed = res_ind_ed(2:2:size(res_ind_ed,1),:);
% % reshape the even rows of the reshaped matrix
% ind_ed2 = [reshape(res_ind_ed',2,[])]';
% %% mod ind_ed with mg to find which micro-grid it belonds
% ind_mod_ed = mod(ind_ed2,mg);
% % if any element is 0 then it is the mg-th micro-grid
% ind_mod_ed(ind_mod_ed == 0) = mg;
% % final matrix of micro-grid scheduling
% mic_mat = ind_mod_ed;
% 
% %% SHORTEST PATHS
% Gs = G0;
% suc_sh_pa = zeros(1,MESS);
% figure(2)
% % h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
% h2 = plot(Gs);
% labelnode(h2,[1:numnodes(Gs)-2],'')
% layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
% title('Shortest Paths')
% disp('******* Shortest Paths Costs ***********')
% for i_rm = 1:MESS
%     [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*','Method','mixed');
% %     [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*');
%     suc_sh_pa(i_rm) = path_len;
%     Rm_nodes = P_nodes(2:(end-1));
% %     G = rmnode(G,Rm_nodes);
%     Gs.Edges.Weight(path1) = inf;
%     highlight(h2,'Edges',path1,'EdgeColor','r','LineWidth',1.5)
%     layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
%     disp(['Cost of ',num2str(i_rm),' path is: ',num2str(suc_sh_pa(i_rm))])
%     %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%     % each column of Shortest_Path_mat is the edges index for each path 
%     Shortest_path_mat(:,i_rm) = path1;
% end
% suc_sh_pa;
% disp(['Total cost reduction with shortest paths is  ',num2str(sum(suc_sh_pa))])
% disp(['Total cost  with shortest paths is  ',num2str(sum(suc_sh_pa)+Base_cost)])
% disp(newline)
% %% DEBUGGING check if the paths are the same (LP & shortest paths)
% disp('**************** Debug Section ************')
% if isequal(fl,isinf(Gs.Edges.Weight))== 1
%     disp('Paths are the same')
% else
%     disp('Paths are not the same')
% end
% if sum(suc_sh_pa) == (G0.Edges.Costs'*fl) 
%     disp('Costs are equal')
% else
%     disp('Costs are not equal')
% end
% disp('******************************')
% disp(newline)
% %% find the paths from LP min cost flow solution 
% % ind_ed_s = ed_mat((fl1>0),:);
% % start and end nodes of edges with flow
% [st_nd, end_nd] = findedge(G0,find(fl1>0));
% ind_ed_s = [st_nd end_nd];
% %  number of edges between S* and T* w/o weights
% dis_S_T = distances(G0,numnodes(G0)-1,numnodes(G0),'Method','unweighted');
%  % distance between day 1 and last day:
%  dis_d1_dl = dis_S_T - 2;
% %%
% [st_nd, end_nd] = findedge(G0,find(fl1>0));
% ind_ed_p = [st_nd end_nd];
% COUNTER = 1;
% init_nodes = ind_ed_p(1:MESS,1);
% %% add the first edges of each path 
% for j_MESS = 1:MESS
% % end_nd_temp = ind_ed_p(COUNTER,2);
% end_nd_temp = ind_ed_p(1,2);
% % end_nd_temp = ind_ed_p(1,1);
% for i_path = 1:(dis_d1_dl-1)
%     node_array(i_path,COUNTER) = end_nd_temp;
% %     outedges(G0,end_nd_temp)
%     next_edge = fl(outedges(G0,end_nd_temp))'*outedges(G0,end_nd_temp);
%     next_edge_array(i_path,COUNTER) = next_edge;
%     % successors(G0,end_nd_temp)
%     [st_nd_temp, end_nd_temp] = findedge(G0,next_edge);
% %     node_mat(i_path,COUNTER) = end_nd_temp;
%     % next node is end_nd_temp
% %     if days <= 7
% %         disp(['Next node is: ',num2str(end_nd_temp)])
% %     end
% end
% node_array(i_path+1,COUNTER) = end_nd_temp;
% next_edge_array_mat(:,j_MESS) = next_edge_array';
% node_array_mat(:,j_MESS) = node_array';
% %
% % find(ind_ed_p(:,1) == node_mat')
% % Give index if the element in "node_mat" is a member of "ind_ed_p(:,1)".
% [~, del_ind2] = ismember(node_array, ind_ed_p(:,2));
% % delete the path from the ind_ed_p matrix and repeat
% ind_ed_p(del_ind2,:) = [];
% % ind_ed_p(del_ind2(2:end),:) = [];
% end
% % node_array_mat(end+1,:) = ind_ed_p(:,2)
% % insert init nodes to the top of the node matrix
% node_array_mat(end+1,:) = init_nodes;
% node_array_mat = circshift(node_array_mat,1);
% % find the edges for the starting nodes
% init_edges = findedge(G0,node_array_mat(1,:),node_array_mat(2,:))';
% % append init edges to the top of next_edge_array_mat
% next_edge_array_mat(end+1,:) = init_edges;
% next_edge_array_mat = circshift(next_edge_array_mat,1);
% %%
% disp(['Cost of each path with LP is',newline, num2str(sum(G0.Edges.Costs(next_edge_array_mat)))])
% path_cost = sum(G0.Edges.Costs(next_edge_array_mat));
% disp(['Total Cost of  LP is',newline, num2str(sum(sum(G0.Edges.Costs(next_edge_array_mat))))])
% %%
% %% DEBUG compare pahts from 2 solutions
% % remove 1st and last row of Shortest path mat cuz it's the edges S* T*
% % if condition holds to run the code more than once
% if size(Shortest_path_mat,1) ~= size(next_edge_array_mat,1)
%     Shortest_path_mat(1,:) = [];
%     Shortest_path_mat(end,:) = [];
% end
% figure(100)
% plot([sum([cumsum(G0.Edges.Costs(next_edge_array_mat))],2) sum(cumsum(Gs.Edges.Costs(Shortest_path_mat)),2)])
% % title('Comparison increasing path cost')
% xlabel('Days')
% ylabel('Cost ($)')
% xticks(1:2:2*days)
% xticklabels(1:days)
% xlim([1 2*days])
% y_min = min(min([sum([cumsum(G0.Edges.Costs(next_edge_array_mat))],2) sum(cumsum(Gs.Edges.Costs(Shortest_path_mat)),2)]));
% ylim([y_min inf])
% legend('LP','Greedy','Location','northwest')
% set(gca,'YGrid','on')