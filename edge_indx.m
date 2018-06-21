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
%%
highlight(h55,'Edges',dly_op_cost_ind(:,day_no),'EdgeColor','r','LineWidth',1.5)
%%
reloc_edges = hglht_edges;
reloc_edges(:,(1:(mg+1):size(hglht_edges,2))) = []
% highlight(h55,'Edges',reloc_edges,'EdgeColor','r','LineWidth',1.5)
%%
highlight(h55,'Edges', reloc_edges(:,4),'Edgecolor','r')
%%
% relocation cost is 0: MESS stays on same micro-grid
highlight(h55,'Edges', (mg+1):(mg+1):max(max(reloc_edges)),'Edgecolor','r')
%%
highlight(h55,'Edges', edge_hlt_tst,'Edgecolor','r')
%%
diag_counter = 1;
day_counter = 5;
% day_counter = max(max(reloc_edges)) - diag_counter;
[diag(reloc_edges(:,1:mg),-diag_counter) ; diag(reloc_edges(:,1:mg),diag_counter)]
dist_mat = [diag(reloc_edges(:,1:mg),-diag_counter) ; diag(reloc_edges(:,1:mg),diag_counter)]+(0:day_counter)*mg*(mg+1)
edge_hlt_tst = reshape(dist_mat,1,[])
% 5+dis_counter*mg
%% naming the nodes to v_11 
rsh_nd_mat  = reshape((1:numnodes(G0)-2),mg,[]);
[I_ndnm,J_ndnm] = find(rsh_nd_mat>0);
str_array = [num2str(round(I_ndnm)) num2str(round(J_ndnm/2))];
num_array = str2num(str_array);
num_mat = reshape(num_array,mg,[]);
num_str = string(num_mat)