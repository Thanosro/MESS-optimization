%% FIGURE 8 MESS-ESS/MESS
mg = 100; days =10;
for MESS = 1:mg
    disp(['------------------sim ',num2str(MESS),'--------------'])
    run ESS_comp.m
    LP_cost(MESS,:) = cost_v*fl ;
    ESS_cost_mat(MESS,:) = sum(ESS_cost);
    disp(newline)
end
save ESS_MESS_ST_kWh.mat mg days  ESS_cost_mat LP_cost
%%
load ESS_MESS_results_kWh.mat
%%
figure(99)
plot([abs(LP_cost)])
hold on
plot(abs(ESS_cost_mat))
xlabel('No of MESS')
ylabel('Gain ($)')
grid on
legend('MESS','ESS','Location','northwest')
%%
close all
%%
% ESS_cost_mat2 = [-6604.6094;-13186.1773
% ;-19764.4718;-26322.272;-32857.5998;-39392.6124;]
figure(1123)
plot([abs(LP_cost)-abs(ESS_cost_mat)])
xlabel('No of MESS')
ylabel('Gain Difference ($)')
grid on
%%
figure(990)
subplot(2,1,1)
plot([abs(LP_cost)])
hold on
plot(abs(ESS_cost_mat))
xlabel('(a)')
ylabel('Gain ($)')
grid on
legend('MESS','ESS','Location','northwest')
subplot(2,1,2)
%
% ESS_cost_mat2 = [-6604.6094;-13186.1773
% ;-19764.4718;-26322.272;-32857.5998;-39392.6124;]
% figure(1123)
plot([abs(LP_cost)-abs(ESS_cost_mat)])
xlabel(['No of MESS',newline,'(b)'])
ylim([0 14000])
ylabel('Gain Difference ($)')
grid on
%%
load ESS_MESS_ST_kWh.mat
figure(11443434)
plot(100*[abs(LP_cost)-abs(ESS_cost_mat)]./abs(ESS_cost_mat))
xlabel('No of MESS')
% title('(MESS-ESS)/ESS')
ylabel('Gain Percentage')
ax = gca;
ax.YGrid = 'on';
ylim([0 35])
% percentage of improvement in benefit using MESS vs ESS decreases 
%%
[max_Gain,max_ind] = max([abs(LP_cost)-abs(ESS_cost_mat)])
%%
print -depsc2 Fig_8_MESS_ESS.eps
%%
save ESS_MESS_results_25.mat mg days max_Gain max_ind ESS_cost_mat LP_cost a_price p_price