clc; clear all;
x = [0 0 1;0 1 0;0 1 0];
Ns = 2; O = 3; T = 4;
x(:,:,1) = x;
x(:,:,2) = [0 1 0;0 1 0;0 0 0];
x(:,:,3) = [1 0 0;0 1 0;1 0 0]
gamma = randi(4,3,3,3)
% gamma = randi(4,Ns,O,T-1)
% electricity price
a = 0.15;
% discharge rate
d = 
% MESS capacity 
s_cap = 
% gamma.^2
% objective cost function
c = a*gamma.^2
c.*x
sum(sum(sum(c.*x)))
% st battery storage 
s_t_mat = randi(4,1,T-1)
s_t = sum(s_t_map)
% constraint #3
% sum of all rows of matrix x
sum_mat = sum(x,2);
if  all(all(sum(x,2) == 1))
    disp('True')
else
    disp("tzifos")
end

% constraint 4
s_t(:) < d
cumsum(s_t) < s_cap