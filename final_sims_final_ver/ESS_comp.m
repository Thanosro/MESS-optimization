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
% G0.Edges.Costs((numedges(G0)-2*mg+1):numedges(G0)) = 0;
G0.Edges.Costs((numedges(G0)-2*mg+1):numedges(G0)) = mean(reloc_mat(:,1));
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