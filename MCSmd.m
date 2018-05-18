% Load_dem = load demand vector
% No_gen_units = number of generation units Ns
% Supply_vec = supply vector of each gen. unit (gamma) (matrix Ns*O*T)
% Ns rows, O columnts T height
% Cost_fun = quadratic function of energy generation Cost_fun =
% a*Supply_vec^2; a = 0.15 $/Kwh
% d upper & lower bound of discharge rate of MESS
% T total time intervals
% x binary vector; shows the operational state of each gen. unit
% O number of operational states of each generation unit
%
function cost_stor = MCSmd(Load_dem,No_gen_units,Supply_vec)
    a = 0.15; % usd/kwh
    Cost_vec = a.*Supply_vec;
    cvx_begin quiet
    cvx_solver gurobi
    variable x(O) binary 
    minimize sum(sum(sum(MATRIX)))
    subject to 
        
    cvx_end
end