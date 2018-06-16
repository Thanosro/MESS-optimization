% find(fl>0)
[st_nd, end_nd] = findedge(G0,find(fl1>0));
ind_ed_p = [st_nd end_nd];
COUNTER = 1;
end_nd_temp = ind_ed_p(COUNTER,2);
%%
for i_path = 1:(dis_d1_dl-1)
    outedges(G0,end_nd_temp)
    next_edge = fl(outedges(G0,end_nd_temp))'*outedges(G0,end_nd_temp);
    next_edge_array(i_path,COUNTER) = next_edge;
    % successors(G0,end_nd_temp)
    [st_nd_temp, end_nd_temp] = findedge(G0,next_edge);
    node_mat(i_path,COUNTER) = end_nd_temp;
    % next node is end_nd_temp
    disp(['Next node is: ',num2str(end_nd_temp)])
end
next_edge_array
node_mat
%%

COUNTER = COUNTER + 1;