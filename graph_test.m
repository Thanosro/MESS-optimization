% block matrix with ones(m)
clc; clear all; close all;
%m # of MESS, d # of days
m=5; d = 7;
% create the adjacency matrix Ad_g
%1st: create block diagonal matrix with ones sub-blocks
Ad_g = kron(eye(d-1),ones(m));
%2nd: include zeros in the 1st m columns of the matrix
Ad_g = [zeros(m*d-m,m) Ad_g];
%3ed: include zeros in the last m rows of the matrix
Ad_g = [Ad_g; zeros(m,m*d)]
% G is the graph
G = digraph(Ad_g);
% plot the graph
plot(G)
