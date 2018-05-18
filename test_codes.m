clc; clear all;
x = [0 0 1;0 1 0;0 1 0];
Ns = 2; O = 3; T = 4;
x(:,:,1) = x;
x(:,:,2) = [0 1 0;0 1 0;0 0 1];
x(:,:,3) = [1 0 0;0 1 0;1 0 0]
a = randi(4,3,3,3)
% a = randi(4,Ns,O,T)

% a.^2
% cost function
c = 0.15*a.^2
c.*x
sum(sum(sum(c.*x)))