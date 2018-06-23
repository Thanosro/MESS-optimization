for i_SP_LP = 1:1000
    run LP_SP_comp.m
    LP_path_cost(i_SP_LP,:) = [path_cost sum(sum(G0.Edges.Costs(next_edge_array_mat)))] ;
    SP_path_cost(i_SP_LP,:) = [suc_sh_pa sum(suc_sh_pa)];
end
%%
PATH_NUMBER = 4;
figure(133)
plot([LP_path_cost(:,PATH_NUMBER) SP_path_cost(:,PATH_NUMBER)])
xlabel('No. of Simulations')
ylabel('Cost ($)')
legend('LP','Greedy','Location','northwest')