% find(fl>0)
[st_nd, end_nd] = findedge(G0,find(fl1>0));
ind_ed_p = [st_nd end_nd];
COUNTER = 1;
init_nodes = ind_ed_p(1:MESS,1);
%% add the first edges of each path 
for j_MESS = 1:MESS
% end_nd_temp = ind_ed_p(COUNTER,2);
end_nd_temp = ind_ed_p(1,2);
% end_nd_temp = ind_ed_p(1,1);
for i_path = 1:(dis_d1_dl-1)
    node_array(i_path,COUNTER) = end_nd_temp;
%     outedges(G0,end_nd_temp)
    next_edge = fl(outedges(G0,end_nd_temp))'*outedges(G0,end_nd_temp);
    next_edge_array(i_path,COUNTER) = next_edge;
    % successors(G0,end_nd_temp)
    [st_nd_temp, end_nd_temp] = findedge(G0,next_edge);
%     node_mat(i_path,COUNTER) = end_nd_temp;
    % next node is end_nd_temp
    disp(['Next node is: ',num2str(end_nd_temp)])
end
node_array(i_path+1,COUNTER) = end_nd_temp;
next_edge_array_mat(:,j_MESS) = next_edge_array';
node_array_mat(:,j_MESS) = node_array';
%
% find(ind_ed_p(:,1) == node_mat')
% Give index if the element in "node_mat" is a member of "ind_ed_p(:,1)".
[~, del_ind2] = ismember(node_array, ind_ed_p(:,2));
% delete the path from the ind_ed_p matrix and repeat
ind_ed_p(del_ind2,:) = [];
% ind_ed_p(del_ind2(2:end),:) = [];
end
% node_array_mat(end+1,:) = ind_ed_p(:,2)
% insert init nodes to the top of the node matrix
node_array_mat(end+1,:) = init_nodes;
node_array_mat = circshift(node_array_mat,1);
% find the edges for the starting nodes
init_edges = findedge(G0,node_array_mat(1,:),node_array_mat(2,:))';
% append init edges to the top of next_edge_array_mat
next_edge_array_mat(end+1,:) = init_edges;
next_edge_array_mat = circshift(next_edge_array_mat,1);
%%
disp(['Cost of each path with LP is',newline, num2str(sum(G0.Edges.Costs(next_edge_array_mat)))])
% sum(G0.Edges.Costs(next_edge_array_mat))
disp(['Total Cost of  LP is',newline, num2str(sum(sum(G0.Edges.Costs(next_edge_array_mat))))])
%%
length(ind_ed_p)
%%

% COUNTER = COUNTER + 1;