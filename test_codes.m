clc; clear all;
% x = [0 0 1;0 1 0;0 1 0];
% Number of generation units
Ns = 2;
% generation operational states
O = 5;
% Time-slots within day
T = 4;
% maximum load demand in any timeslot (kW)
MAX_DEM = 100;
% maximum energy gen. in any timeslot (kW)
MAX_GEN = 100;
% electricity price (usd $)
a = 0.15;
% discharge rate (kW)
d = 35;  %the mess discharges 4 kW per hour
% MESS max capacity (kW)
s_cap = 40; 
x_2d = [0 1 0 0 0;0 0 0 1 0];
x_3d = repmat(x_2d,1,1,T-1);

% power generation (kW) 
% gamma_1d = [0 25 50 100]
% gamma_1d = 0:25:100
gamma_1d = linspace(0,MAX_GEN,O); % distinct energy gen. values (0 25 50 100)
gamma_2d = repmat(gamma_1d,Ns,1); % number of generators arrays
gamma_3d = repmat(gamma_1d,Ns,1,T-1); % Time slots arrays

% Load demand L(t) (kW)
% array with the load demand each time-slot  
L_t_ar =  randi(MAX_DEM,1,T-1); 
% gamma.^2
% objective cost function
c = a*gamma_3d.^2;
c.*x_3d
obj_fun = sum(sum(sum(c.*x_3d)))
% st battery storage 
s_t_mat = randi(d,1,1,T-1);
s_t = sum(s_t_mat);
% constraint #2
L_t_ar <= sum(sum(gamma_3d(:,:,:).*x_3d(:,:,:))) + s_t_mat
% constraint #3
% sum of all rows of matrix x
sum_mat = sum(x_3d,2);
if  all(all(sum(x_3d,2) == 1))
    disp('True')
else
    disp("tzifos")
end
% constraint 4
s_t_mat(:) <= d
cumsum(s_t_mat(:)) <= s_cap
%%
% clc; clear all;
% Number of generation units
Ns = 2;
% generation operational states
O = 5;
% Time-slots within day
T = 4;
% maximum load demand in any timeslot (kW)
MAX_DEM = 100;
% maximum energy gen. in any timeslot (kW)
MAX_GEN = 50;
% L_t = randi(MAX_DEM,1,1,T-1);
L_t = MAX_DEM*ones(1,1,T-1);
gamma_1d = linspace(0,MAX_GEN,O); % distinct energy gen. values (0 25 50 100)
% gamma_2d = repmat(gamma_1d,Ns,1); % number of generators arrays
gamma_3d = repmat(gamma_1d,Ns,1,T-1); % Time slots arrays
[min_cost, stor, x_var] = MCSmd(L_t,Ns,O,T,gamma_3d);
%%
% code eval
% generation units operational strategies
clc;
disp(['Operational Strategies used:'])
op_strategies
disp('Cost of generated power:')
cost_en_3d = 0.25*op_strategies.^2
disp(['Load demand each T is:' newline num2str(L_t)])
disp(['Total cost of gen. power each T is:' newline num2str(sum(sum((cost_en_3d))))])
disp(['Total cost of gen. power is:' newline num2str(sum(sum(sum((cost_en_3d)))))])
disp(['Storage unit used:', newline num2str(stor)])
disp(['Power cost of MESS is:' newline num2str(0.05*((sum(stor))^2)) ])
disp(['Min cost is ',newline num2str(min_cost)])
disp(['Charge left is:',newline num2str(50-sum(stor))])
