clc; clear all;
x = [0 0 1;0 1 0;0 1 0];
Ns = 2; O = 3; T = 4;
x(:,:,1) = x;
x(:,:,2) = [0 1 0;0 1 0;0 0 0];
x(:,:,3) = [1 0 0;0 1 0;1 0 0]
a = randi(4,3,3,3)
% a = randi(4,Ns,O,T-1)

% a.^2
% objective cost function
c = 0.15*a.^2
c.*x
sum(sum(sum(c.*x)))
% constraint #3
% sum of all rows of matrix x
sum_mat = sum(x,2);
if  all(all(sum(x,2) == 1))
    disp('True')
else
    disp("tzifos")
end
