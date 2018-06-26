for i_SP_LP = 1:10%50
    disp(['------------------sim ',num2str(i_SP_LP),'--------------'])
    run LP_SP_comp.m
    LP_path_cost(i_SP_LP,:) = [path_cost sum(sum(G0.Edges.Costs(next_edge_array_mat)))] ;
    SP_path_cost(i_SP_LP,:) = [suc_sh_pa sum(suc_sh_pa)];
    base_cost_mat(i_SP_LP) = Base_cost;
end
%%
% load 30_sims.mat
%%
figure(500)
PATH_NUMBER = 4;
% subplot(2,2,PATH_NUMBER)
plot([LP_path_cost(:,PATH_NUMBER)+base_cost_mat' SP_path_cost(:,PATH_NUMBER)+base_cost_mat' base_cost_mat'])
xlabel(['No. of Simulations'])%,newline,'(',num2str(PATH_NUMBER),')'])
ylabel('Cost ($)')
if PATH_NUMBER == MESS+1
    title('Total Cost')
else
title(['MESS # ',num2str(PATH_NUMBER)])
end
legend('LP','Greedy','Location','northwest')
%%
figure(501)
PATH_NUMBER = 4;
% subplot(2,2,PATH_NUMBER)
plot([LP_path_cost(:,PATH_NUMBER)-SP_path_cost(:,PATH_NUMBER)])
xlabel(['No. of Simulations'])%,newline,'(',num2str(PATH_NUMBER),')'])
ylabel('Cost ($)')
if PATH_NUMBER == MESS+1
    title('Total Cost')
else
title(['MESS # ',num2str(PATH_NUMBER)])
end
legend('LP','Greedy','Location','northwest')

%%
% for PATH_NUMBER = 1:MESS+1
    figure(445)
    subplot(2,2,PATH_NUMBER)
    plot([LP_path_cost(:,PATH_NUMBER) SP_path_cost(:,PATH_NUMBER)])
    xlabel(['No. of Simulations',newline,'(',num2str(PATH_NUMBER),')'])
    ylabel('Cost ($)')
    if PATH_NUMBER == MESS+1
        title('Total Cost')
    else
    title(['MESS # ',num2str(PATH_NUMBER)])
    end
    legend('LP','Greedy','Location','northwest')
% end
%%

figure(1344)
plot([-LP_path_cost(:,4)+SP_path_cost(:,4)])
xlabel('No. of Simulations')
ylabel('Gain ($)')
% legend('Gain','Location','northwest')
%%
SIM_NO = 1;
figure(166)
plot([LP_path_cost(SIM_NO,1:3)])
hold on
plot([SP_path_cost(SIM_NO,1:3)])%; SP_path_cost(SIM_NO,1:3)])
plot([LP_path_cost(SIM_NO,1:3) - SP_path_cost(SIM_NO,1:3)])
xlabel('No. of Simulations')
ylabel('Gain ($)')


%%
% save 34_sims.mat LP_path_cost SP_path_cost mg days MESS G0
LP_mean = mean(LP_path_cost)
SP_mean = mean(SP_path_cost)
gain = SP_mean - LP_mean
gain_per = gain(4)/LP_mean(4)