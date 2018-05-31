% graph with node split
clc; clear all; close all;
%m , d # of days
mg = 3; days = 3;
% # of MESS
K = 2;
if K>mg
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
G0 = addnode(G0,{'S*'});
G0 = addnode(G0,{'T*'});
Ed_S0 = table([(2*mg*days+1)*ones(mg,1) (1:mg)'],K*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T0 = table([ ((2*mg*days-mg+1):(2*mg*days))' (2*mg*days+2)*ones(mg,1)],K*ones(mg,1),'VariableNames',{'EndNodes','Weight'});
G0 = addedge(G0,Ed_S0);
G0 = addedge(G0,Ed_T0);
% G0 = addedge(G0,2*mg*days+2,2*mg*days+1,K);
figure(232)
h3 =plot(G0,'EdgeLabel',G0.Edges.Weight);
layout(h3,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% highlight(h3,'Edges',,'EdgeColor','r')
%%
cost_v = randi(5,1,numedges(G0));
%%
cvx_begin
cvx_solver gurobi
variable fl(numedges(G0),1)
% maximize sum(cost_v*fl)
minimize cost_v*fl
subject to 
%     incidence(G0)*fl == 0;
    incidence(G0)*fl == [zeros(numnodes(G0)-2,1); -K; K];
    0 <= fl <= G0.Edges.Weight;
cvx_end
[fl cost_v']
%%
highlight(h3,'Edges',find(fl>0),'EdgeColor','r','LineWidth',1.5)