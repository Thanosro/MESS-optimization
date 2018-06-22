% clc; clear all;
% load the array 
% loads as struct
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
tic
% for i_mcs = 1:mg*days
i_mcs = 1;
while i_mcs <=  mg*days
    disp(num2str(i_mcs))
cl_num = randi(9000);
L_t_96 = 25*L_t_ar(96*cl_num:96*(cl_num+1)-1);
if any(isnan(L_t_96))>0
    cl_num = randi(9000);
    L_t_96 = 25*L_t_ar(96*cl_num:96*(cl_num+1)-1);
    disp('Nan')
end
% L_t_96 has 96 values, each 15 mins
% convert to hourly load demands
L_t_mat = reshape(L_t_96,4,T);
L_t_h = sum(L_t_mat);
L_t_h(T) = [];
L_t_h3d = reshape(L_t_h,1,1,T-1);
[min_cost_S, stor_S, x_var_S] = MCSmd(L_t_h3d,Ns,O,T,gamma_3d);
[min_cost, stor, x_var] = MCmd(L_t_h3d,Ns,O,T,gamma_3d);
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
%
d_cost_red = diff(min_cost_out,1,2);
% reshape(d_cost_red,mg,[])
% d_cost_red = reshape(dly_op_cost_ind,[],1);
%%
cl_num = randi(9000);
L_t_96 = L_t_ar(96*cl_num:96*(cl_num+1)-1);
L_t_mat = reshape(L_t_96,4,T);
L_t_h = sum(L_t_mat);
L_t_h(T) = [];
L_t_h3d = reshape(L_t_h,1,1,T-1);
tot_dail_load = sum(sum(sum(L_t_96)))

