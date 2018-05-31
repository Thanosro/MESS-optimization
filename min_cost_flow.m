% syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 real 
% x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]';
x = randi(4,numedges(G),1);
incidence(G)*x
%%
cost_v = randi(5,1,numedges(G));
%%
cvx_begin
cvx_solver gurobi
variable fl(numedges(G),1)

% maximize sum(cost_v*fl)
maximize cost_v*fl
subject to 
    incidence(G)*fl == 0;
    0 <= fl <= G.Edges.Weight;
cvx_end
[fl cost_v']
 %% show edges on graph plot
highlight(h1,'Edges',find(fl>0),'EdgeColor','r','LineWidth',1.5)