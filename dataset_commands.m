addpath(genpath('C:\Users\Thanos\Documents\DeepSolar'))
cd C:\Users\Thanos\Documents\DeepSolar\Systech\sims\MESS-optimization
%% read dataset commands
L_t_ar = xlsread('gen_2015.csv','C:C');
%%

plot(L_t_ar(96*cl_num:96*(cl_num+1)))