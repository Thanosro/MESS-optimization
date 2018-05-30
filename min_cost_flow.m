syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 real 
x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13]';
incidence(G)*x
%%
cost_v = randi(5,1,numedges(G));
cvx_begin
cvx_solver gurobi
variable fl(numedges(G))
minimize cost_v*fl
subject to 
    incidence(G)*fl == 0;
   0<= fl <= G.Edges.Weight;
cvx_end
%% show edges on graph plot
highlight(h,find(fl>0),'EdgeColor','r')
