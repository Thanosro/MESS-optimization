% FIGURE 6  paths for K=1,2,4
addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
rmpath('C:\Users\Thanos\Documents\DeepSolar\Optimal_flow\cvx\lib\narginchk_')
%% Laptop
addpath(genpath('C:\Users\thano\OneDrive\Documents\USC'))
cd C:\Users\thano\OneDrive\Documents\USC\DeepSolar\sustech\MESS-optimization-master\MESS-optimization
rmpath('C:\Users\thano\OneDrive\Documents\USC\DeepSolar\OPF\cvx\lib\narginchk_')
%% close all;
%m , d # of days
mg = 4; days =4;
% # of MESS
%%
L_t_str = load('Pecan_load.mat');
% convert to array 
L_t_ar = L_t_str.L_t_ar;
% Number of generation units
Ns = 4;
% generation operational states
O = 5;
% Time-slots within day
T = 24;
% mg = 3; days = 4;
% maximum load demand in any timeslot (kW)
MAX_DEM = 100;
% maximum energy gen. in any timeslot (kW)
MAX_GEN = 150;
% L_t = randi(MAX_DEM,1,1,T-1);
% L_t = MAX_DEM*ones(1,1,T-1);
gamma_1d = linspace(0,MAX_GEN,O); % distinct energy gen. values (0 25 50 100)
% gamma_2d = repmat(gamma_1d,Ns,1); % number of generators arrays
gamma_3d = repmat(gamma_1d,Ns,1,T-1); % Time slots arrays
%%
stor_cap = 120; % kWh
dis_rate = 30;
%%
tic
% for i_mcs = 1:mg*days
i_mcs = 1;
while i_mcs <=  mg*days
    disp(num2str(i_mcs))
cl_num = randi(5000);
L_t_96 =  25*L_t_ar(96*cl_num:96*(cl_num+1)-1);
if any(isnan(L_t_96))>0
    cl_num = randi(5000);
    L_t_96 = 25*L_t_ar(96*cl_num:96*(cl_num+1)-1);
    disp('Nan')
end
% L_t_96 has 96 values, each 15 mins
% convert to hourly load demands
L_t_mat = reshape(L_t_96,4,T);
L_t_h = mean(L_t_mat);
L_t_h(T) = [];
L_t_h3d = reshape(L_t_h,1,1,T-1);
[min_cost_S, stor_S, x_var_S] = MCSmd(L_t_h3d,Ns,O,T,gamma_3d,stor_cap,dis_rate);
[min_cost, stor, x_var] = MCmd(L_t_h3d,Ns,O,T,gamma_3d,stor_cap);
min_cost_out(i_mcs,:) = [min_cost_S min_cost];
if any(any(isnan(min_cost_out(i_mcs,:))))>0
    i_mcs = i_mcs - 1;
    disp('Infeasible')
end
i_mcs = i_mcs + 1;
end
tim = toc
disp('done')
min_cost_out;
% cost of nodes
disp(['Base operational cost without MESS is: ',num2str(sum(min_cost_out(:,2)))])
Base_cost =  sum(min_cost_out(:,2));
d_cost_red = -diff(min_cost_out,1,2);
%%
MESS_MAT = [1 2 4];
mess_counter = 3;
MESS = MESS_MAT(mess_counter);
if MESS > mg
       error('No of MESS > of Micro-grids')
end
% eye(mg) is the diagonal with the cost difference with & w/o MESS
% ones is the matrix denoting the tranfer cost from microgrid i to j
% create the relocation cost matrix
base_reloc_cost = 5;
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
figure(1001)
subplot(3,2,2*mess_counter-1)
% h1 =plot(G0,'EdgeLabel',G0.Edges.Costs);
h1 =plot(G0);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
if MESS == 1
    title('Min cost flow LP','FontSize',12,'FontWeight','bold')
end
labelnode(h1,[1:numnodes(G0)-2],'')
ylabel(['K = ',num2str(MESS)],'FontSize',12,'FontWeight','bold')
xlabel('Days')%,'FontSize',10,'FontWeight','bold')
set(gca,'YTick',[])
set(gca,'XTick',[])
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
disp(['Total cost reduction with LP: ',num2str(cost_v*fl)])
disp(['Total cost   with LP: ',num2str(cost_v*fl+Base_cost)])
%%
highlight(h1,'Edges',find( fl>0),'EdgeColor','r','LineWidth',1.5)
%% SHORTEST PATHS
Gs = G0;
suc_sh_pa = zeros(1,MESS);
% figure(2)
subplot(3,2,2*mess_counter)
% h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
h2 = plot(Gs);
labelnode(h2,[1:numnodes(Gs)-2],'')
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
if MESS == 1
    title('Shortest Paths','FontSize',12,'FontWeight','bold')
end
xlabel('Days')%,'FontSize',10,'FontWeight','bold')
disp('******* Shortest Paths Costs ***********')
for i_rm = 1:MESS
    [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*','Method','mixed');
%     [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*');
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
disp(['Total cost reduction with shortest paths is  ',num2str(sum(suc_sh_pa))])
disp(['Total cost  with shortest paths is  ',num2str(sum(suc_sh_pa)+Base_cost)])
disp(newline)
set(gca,'YTick',[])
set(gca,'XTick',[])
%%
print -depsc2 myplot.eps

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