% syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 real 
% x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]';
x = randi(4,numedges(G),1);
incidence(G)*x
%%
cost_v = randi(5,1,numedges(G));
%%
cost_vec = G.Edges.Costs';
cvx_begin
cvx_solver gurobi
variable fl(numedges(G),1)

% maximize sum(cost_v*fl)
<<<<<<< Updated upstream
minimize cost_vec*fl
subject to 
    incidence(G)*fl == [ zeros((numnodes(G)-2),1); -MESS; MESS];
=======
maximize cost_vec*fl
subject to 
    incidence(G)*fl == zeros(14,1);
>>>>>>> Stashed changes
    0 <= fl <= G.Edges.Capacities;
cvx_end
[fl cost_vec']
 %% show edges on graph plot
highlight(h1,'Edges',find(fl>0),'EdgeColor','r','LineWidth',1.5)
%% find path of nodes 
G.Edges.EndNodes([fl>0 fl>0])