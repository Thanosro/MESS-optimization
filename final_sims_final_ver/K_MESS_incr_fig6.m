%% greedy vs LP relocation cost figure
% graph with node split
clc; clear all;
clc;
% clearvars -except i_SP_LP LP_path_cost SP_path_cost MESS mg
mg = 100; days =10;
load 100_mg_10d_Kwh.mat
%%
close all;
reloc_factor = 2;
disp(['Cost without MESS is :',num2str(sum(min_cost_out(:,2)))])
% sum(min_cost_out(:,2))
disp(['Cost with MESS in every micro-grid is:',num2str(sum(min_cost_out(:,1)))])
disp(['Reloc factor is ',num2str(reloc_factor)])
for MESS = 1:mg
%m , d # of days
mg = 100; days =10;
% # of MESS
% MESS = 45;
% cost without MESS
disp(['--------------',num2str(MESS),'----------------'])
if MESS > mg
       error('No of MESS > of Micro-grids')
end
% eye(mg) is the diagonal with the cost difference with & w/o MESS
% ones is the matrix denoting the tranfer cost from microgrid i to j
% create the relocation cost matrix
base_reloc_cost = reloc_factor*5;
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
% add edges costs and capacities 
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
% figure(1)
% % h1 =plot(G0,'EdgeLabel',G0.Edges.Costs);
% h1 =plot(G0);
% layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% title('Min cost flow LP')
% labelnode(h1,[1:numnodes(G0)-2],'')
% % highlight(h3,'Edges',,'EdgeColor','r')
%
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
disp(['Total cost with LP  ',num2str(MESS), ' MESS: ',num2str(cost_v*fl)])
disp(newline)
cost_red_mat_LP(MESS) = cost_v*fl;
%######################################
suc_sh_pa = zeros(1,MESS);
Gs = G0;
title('Shortest Paths')
disp('******* Shortest Paths Costs ***********')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*','Method','mixed');
    suc_sh_pa(i_rm) = path_len;
    Rm_nodes = P_nodes(2:(end-1));
    Gs.Edges.Weight(path1) = inf;
%     disp(['Cost of ',num2str(i_rm),' path is: ',num2str(suc_sh_pa(i_rm))])
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % each column of Shortest_Path_mat is the edges index for each path 
    Shortest_path_mat(:,i_rm) = path1;
end
% disp(['Total cost reduction with shortest paths is  ',num2str(sum(suc_sh_pa))])
disp(['Total cost  with shortest paths is  ',num2str(sum(suc_sh_pa))])
disp(newline)
cost_red_mat_SP(MESS) = sum(suc_sh_pa);
clearvars -except cost_red_mat_SP MESS cost_red_mat_LP mg days d_cost_red reloc_factor
end
%%
save K_incr_ST_6_KwH_rlfc_2.mat
%%
% save K_incr_ST_6_KwH_rlfc_1.mat
%%
% save K_incr_fig_6_KwH_rlfc_3.mat
% %%
% save K_incr_fig_6_KwH_rlfc_10.mat
%%
load  K_incr_ST_6_KwH_rlfc_2.mat % reloc factor = 1
%%
load K_incr_ST_6_KwH_rlfc_1.mat % reloc_factor = 2
%%
load  K_incr_ST_6_KwH_rlfc_1.mat % reloc factor = 1
Gain_LP_A1 = abs(cost_red_mat_LP);
Gain_SP_G1 = abs(cost_red_mat_SP);
Gain_perc1 = (Gain_LP_A1 - Gain_SP_G1)./ Gain_SP_G1;
%%
load K_incr_ST_6_KwH_rlfc_2.mat % reloc_factor = 2
Gain_LP_A2 = abs(cost_red_mat_LP);
Gain_SP_G2 = abs(cost_red_mat_SP);
Gain_perc2 = (Gain_LP_A2 - Gain_SP_G2)./ Gain_SP_G2;
%%
figure(10001)
plot(100*Gain_perc1)
% title('(A-G)/G')
hold on
plot(100*Gain_perc2)
xlabel('No of MESS')
ylabel('Gain Percentage')
grid on
ylim([0 inf])
ax = gca;
ax.XGrid = 'off';
ax.YGrid = 'on';
yticks([0 2.5 5 7.5 10 12.5 15])
legend(['R_L = ',num2str(1)],['R_L = ',num2str(2)],'location','northwest')
%%
print -depsc2 Fig6_R_L_ST.eps
%%
figure(33333)
plot([abs(cost_red_mat_LP) ])
hold on
plot([abs(cost_red_mat_SP) ])
% plot(cost_red_mat_SP)
hold on
% plot(suc_sh_pa)
grid on
xlabel('No of MESS')
legend('Min-Cost Flow','Greedy','Location','northwest')
ylabel('Gain ($)')

%%
figure(33366)
subplot(2,1,1)
plot([abs(cost_red_mat_LP) ])
hold on
plot([abs(cost_red_mat_SP) ])
xlabel('(a)')
ylabel('Gain ($)')
grid on
legend('Min-Cost Flow','Greedy','Location','northwest')
subplot(2,1,2)
%
% ESS_cost_mat2 = [-6604.6094;-13186.1773
% ;-19764.4718;-26322.272;-32857.5998;-39392.6124;]
% figure(1123)
plot([abs(cost_red_mat_LP) - abs(cost_red_mat_SP)])
xlabel(['No of MESS',newline,'(b)'])
% ylim([0 3500])
ylabel('Gain Difference ($)')
grid on
iax = 1; % Or whichever
% subaxis(4, 6, iax, 'sh', 0.03, 'sv', 0.01, 'padding', 0, 'margin', 0);
%%
% save K_incr_section_2.mat
%%
% Base_cost_array = sum(reshape(min_cost_out(:,2),mg,[]))
%%
reduced_cost_LP = -abs(cost_red_mat_LP)+sum(min_cost_out(:,2))*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2));
reduced_cost_SP = -abs(cost_red_mat_SP)+sum(min_cost_out(:,2))*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2));
%%
figure(1789)
% plot([abs(cost_red_mat_LP) ])
% hold on
% plot([abs(cost_red_mat_SP) ])

plot(sum(min_cost_out(:,1)*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2))))
hold on
plot(sum(min_cost_out(:,2)*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2))))
hold on
plot(reduced_cost_LP)
hold on
plot(reduced_cost_SP )







