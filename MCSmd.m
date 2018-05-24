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
function cost_stor = MCSmd(Load_dem,No_gen_units,Op_States,Time_slots,Supply_vec)
%     a = 0.15; % usd/kwh
%     Cost_vec = a.*Supply_vec;
% No_gen_units = Ns
% Time_slots = T
% Op_States = O
% Supply_vec = gamma_3d
    cvx_begin quiet
    cvx_solver gurobi
    variable x_3d(No_gen_units,Op_States,Time_slots-1) binary 
    variable s_t_mat(1,1,Time_slots-1)
    minimize sum(sum(sum(MATRIX)))
    subject to 
        % const. #2
        Load_dem = sum(sum(Supply_vec(:,:,:).*x_3d(:,:,:))) + 
        % const. #3
        all(all(sum(x_3d,2))) == 1
        % const. #4
        0<= cumsum(s_t_mat) <= s_cap
        % const. #5 
        s_t_mat(:) <= d
    cvx_end
end