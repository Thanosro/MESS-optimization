addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
%% read dataset commands
% read the 3rd coluumn as an array 
% L_t_ar = xlsread('gen_2015.csv','C:C');
%% load the array 
% loads as struct
L_t_str = load('Pecan_load.mat');
% convert to array 
L_t_ar = L_t_str.L_t_ar;
%% plot clusters of 96 data poits (1 day)
cl_num = 40;
plot(L_t_ar(96*cl_num:96*(cl_num+1)-1))
%%
for i_plot = 1:10
   cl_num = randi(10000);
   96*cl_num;
   figure(131)
   plot(L_t_ar(96*cl_num:96*(cl_num+1)-1))
   hold on;
end
%%
close all;