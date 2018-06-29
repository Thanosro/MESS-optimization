%% figure for percentage gain for K increasin Section V.1 figure 5
load 100_mg_10d_Kwh.mat
load  K_incr_fig_6_KwH.mat
%%
reduced_cost_LP = -abs(cost_red_mat_LP)+sum(min_cost_out(:,2))*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2));
reduced_cost_SP = -abs(cost_red_mat_SP)+sum(min_cost_out(:,2))*ones(size(cost_red_mat_LP,1),size(cost_red_mat_LP,2));

perc_LP = reduced_cost_LP./sum(min_cost_out(:,2));
perc_SP = reduced_cost_SP./sum(min_cost_out(:,2));
%%
figure(11880)
title('Percentage of cost reduction')
plot(100*perc_LP)
xlabel('No of MESS')
ylabel('Total Cost Percentage')
ax = gca;
ax.YGrid = 'on';
%%
print -depsc2 Fig_Kincr_V1.eps