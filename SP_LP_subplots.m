%% MESS = 1
cost_v = G0.Edges.Costs';
tic;
cvx_begin quiet
cvx_solver gurobi
variable fl_11(numedges(G0),1)  %integer
% maximize sum(cost_v*fl)
minimize cost_v*fl_11
subject to 
    incidence(G0)*fl_11 == [zeros(numnodes(G0)-2,1); -MESS; MESS];
    0 <= fl_11 <= G0.Edges.Capacities;
cvx_end
toc;
% [fl cost_v'];
disp(['Total cost reduction with LP: ',num2str(cost_v*fl_11)])
disp(['Total cost   with LP: ',num2str(cost_v*fl_11+Base_cost)])
%%
highlight(h1,'Edges',find( fl_11>0),'EdgeColor','r','LineWidth',1.5)
%%
%% SHORTEST PATHS
Gs = G0;
suc_sh_pa = zeros(1,MESS);
figure(2)
% h2 = plot(Gs,'EdgeLabel',Gs.Edges.Weight);
h2 = plot(Gs);
labelnode(h2,[1:numnodes(Gs)-2],'')
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
title('Shortest Paths')
disp('******* Shortest Paths Costs ***********')
for i_rm = 1:MESS
    [P_nodes_11,path_len_11,path1_11] = shortestpath(Gs,'S*','T*','Method','mixed');
%     [P_nodes,path_len,path1] = shortestpath(Gs,'S*','T*');
    suc_sh_pa(i_rm) = path_len_11;
    Rm_nodes = P_nodes_11(2:(end-1));
%     G = rmnode(G,Rm_nodes);
    Gs.Edges.Weight(path1_11) = inf;
    highlight(h2,'Edges',path1_11,'EdgeColor','r','LineWidth',1.5)
    layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
    disp(['Cost of ',num2str(i_rm),' path is: ',num2str(suc_sh_pa(i_rm))])
    %&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % each column of Shortest_Path_mat is the edges index for each path 
    Shortest_path_mat_11(:,i_rm) = path1_11;
end
suc_sh_pa;
disp(['Total cost reduction with shortest paths is  ',num2str(sum(suc_sh_pa))])
disp(['Total cost  with shortest paths is  ',num2str(sum(suc_sh_pa)+Base_cost)])
disp(newline)
%%
figure(4444)
subplot(2,2,1)
h1 =plot(G0);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
title(['Min-Cost Flow',newline,'K = 2'])
labelnode(h1,[1:numnodes(G0)-2],'')
ylabel('Micro-grids')
xlabel(['Days',newline,'(a)'])
highlight(h1,'Edges',find( fl>0),'EdgeColor','r','LineWidth',1.5)
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,2,2)
h2 = plot(Gs);
labelnode(h2,[1:numnodes(Gs)-2],'')
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
title(['Greedy',newline,'K = 2'])
highlight(h2,'Edges',reshape(Shortest_path_mat,1,[]),'EdgeColor','r','LineWidth',1.5)
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
ylabel('Micro-grids')
xlabel(['Days',newline,'(b)'])
set(gca,'xtick',[])
set(gca,'ytick',[])
% hold on
%
% figure(4444)
subplot(2,2,3)
h1 =plot(G0);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
% title('Min Cost Flow')
labelnode(h1,[1:numnodes(G0)-2],'')
title('K = 1')
ylabel('Micro-grids')
xlabel(['Days',newline,'(c)'])
highlight(h1,'Edges',find( fl_11>0),'EdgeColor','r','LineWidth',1.5)
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,2,4)
h2 = plot(Gs);
labelnode(h2,[1:numnodes(Gs)-2],'')
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
% title('Greedy')
highlight(h2,'Edges',reshape(Shortest_path_mat_11,1,[]),'EdgeColor','r','LineWidth',1.5)
layout(h2,'layered','Direction','right','Sources','S*','Sinks','T*')
title('K = 1')
ylabel('Micro-grids')
xlabel(['Days',newline,'(d)'])
set(gca,'xtick',[])
set(gca,'ytick',[])