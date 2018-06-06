addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
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
G0.Edges.Weight = G0.Edges.Costs;
% assign label nodes
G0.Edges.Labels = (1:numedges(G0))';
figure(1)
h1 =plot(G0,'EdgeLabel',G0.Edges.Costs);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% highlight(h3,'Edges',,'EdgeColor','r')
%%
cost_v = G0.Edges.Costs';
cvx_begin
cvx_solver gurobi
variable fl(numedges(G0),1)
% maximize sum(cost_v*fl)
minimize cost_v*fl
subject to 
    incidence(G0)*fl == [zeros(numnodes(G0)-2,1); -MESS; MESS];
    0 <= fl <= G0.Edges.Capacities;
cvx_end
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


