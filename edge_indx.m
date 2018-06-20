G0.Edges
%%
figure(8)
h55 =plot(G0,'EdgeLabel',G0.Edges.Labels);
layout(h55,'layered','Direction','right','Sources', 'S*','Sinks','T*')
title('LP indices')
% highlight(h55,'Edges',1:mg,'EdgeColor','r','LineWidth',1.5)
%%
% numedges(G0)-2*mg
% it = 2
for it = 0:days-1
    disp(['count: ',num2str(it)])
(it*(mg+(mg^2))+1):(it*(mg+(mg^2))+mg)
highlight(h55,'Edges',(it*(mg+(mg^2))+1):(it*(mg+(mg^2))+mg),'EdgeColor','r','LineWidth',1.5)
end
%% 
hglht_edges = reshape(G0.Edges.Labels,mg,[]);
hglht_edges(:,(end-1:end)) =[];
% each column is the index daily cost of the matrix
dly_op_cost_ind = hglht_edges(:,(1:(mg+1):size(hglht_edges,2)))
%%
day_no = 3 % index of each day 
disp(['The cost of each micro-grid in day ',num2str(day_no),' is:'])
G0.Edges.Costs(dly_op_cost_ind(:,day_no))