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
function [min_cost_S, s_t_mat_3d, x_3d] = MCSmd(Load_dem,No_gen_units,Op_States,Time_slots,Supply_vec_3d)
    a = 0.25; % usd/kwh
    p = 0.05; % usc/kWh (night price)
    d = 20;
    s_cap = 50;
%     Cost_vec = a.*Supply_vec;
% No_gen_units = Ns
% Time_slots = T
% Op_States = O
% Supply_vec = gamma_3d
cost_en_3d = a*Supply_vec_3d.^2;
cvx_begin quiet
cvx_solver gurobi
variable x_3d(No_gen_units,Op_States,Time_slots-1) binary 
variable s_t_mat_3d(1,1,Time_slots-1)
minimize sum(sum(sum(cost_en_3d.*x_3d))) + p*(sum(s_t_mat_3d)^2)
subject to 
    % const. #2
    Load_dem <= sum(sum(Supply_vec_3d(:,:,:).*x_3d(:,:,:))) + s_t_mat_3d;
    % const. #3
    sum(x_3d,2) == 1;
    % const. #4
    0<= cumsum(s_t_mat_3d) <= s_cap;
    % const. #5 
   0 <= s_t_mat_3d(:) <= d;
cvx_end
min_cost_S = sum(sum(sum(cost_en_3d.*x_3d))) + p*(sum(s_t_mat_3d)^2);
end