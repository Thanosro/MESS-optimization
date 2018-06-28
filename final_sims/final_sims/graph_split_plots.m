% plot graphs for figure 3 in paper
% graph with node split
clc; clear all; close all;
%m , d # of days
mg = 3; days =4;
% # of MESS
base_reloc_cost = 5;
reloc_mat =tril(ones(mg));

reloc_mat = tril(reloc_mat);
reloc_mat =  base_reloc_cost*(reloc_mat+triu(reloc_mat',1));
A0 = [eye(mg) zeros(mg);zeros(mg) reloc_mat];
AD0 = kron(eye(days-1),A0);
% spy(AD0)
AD1 = blkdiag(AD0,eye(mg));
AD2 = [zeros((days-1)*2*mg+mg,mg) AD1; zeros(mg) zeros(mg,(days-1)*2*mg+mg)];
% spy(AD2)
G0 = digraph(AD2);
% initial number of edges
init_edges = numedges(G0);
G0 = addnode(G0,{'S*'});
G0 = addnode(G0,{'T*'});
Ed_S0 = table([(2*mg*days+1)*ones(mg,1) (1:mg)'],ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T0 = table([ ((2*mg*days-mg+1):(2*mg*days))' (2*mg*days+2)*ones(mg,1)],ones(mg,1),'VariableNames',{'EndNodes','Weight'});
G0 = addedge(G0,Ed_S0);
G0 = addedge(G0,Ed_T0);
%% add edges costs and capacities 
% edge capacities are 1 except those of S* and T* with capacity equal to
% the # of MESS
% assign label nodes
G0.Edges.Labels = (1:numedges(G0))';
G0.Edges.Capacities = G0.Edges.Labels;

% figure(1)
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,2)
h1 =plot(G0);
layout(h1,'layered','Direction','right','Sources', 'S*','Sinks','T*')
highlight(h1,'Edges',1:numedges(G0),'LineWidth',1.5)
highlight(h1,1:numnodes(G0),'MarkerSize',4)
xd = get(h1, 'XData');
yd = get(h1, 'YData');


set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'Visible','off')

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

box off
xlabel('(b)')
labelnode(h1,[1:numnodes(G0)],'')
highlight(h1,[1:numnodes(G0)]' ,'NodeColor','r','MarkerSize',6)
% text(h1.XData(end-1)+0.1,h1.YData(end-1),'S*','fontsize',10,'HorizontalAlignment','right')
% text(h1.XData(end)+0.1,h1.YData(end),'T*','fontsize',10)
% % highlight(h1,'Edges',path1,'EdgeColor','r','LineWidth',1.5)

%%
hglht_edges = reshape(G0.Edges.Labels,mg,[]);
hglht_edges(:,(end-1:end)) =[];
% each column is the index daily cost of the matrix
dly_op_cost_ind = hglht_edges(:,(1:(mg+1):size(hglht_edges,2)));
dly_op_cost_ind = reshape(dly_op_cost_ind,[],1);
highlight(h1,'Edges',dly_op_cost_ind,'EdgeColor','r','LineWidth',2.5)
%%
% mg = 2; days =3;
% % # of MESS
% base_reloc_cost = 5;
% reloc_mat =tril(ones(mg));
% 
% reloc_mat = tril(reloc_mat);
% reloc_mat =  base_reloc_cost*(reloc_mat+triu(reloc_mat',1));
A0o = [ones(mg) ];
AD0o = kron(eye(days-1),A0o);
% spy(AD0)
% AD1o = blkdiag(AD0o,eye(mg));
AD1o = [zeros((days-1)*mg,mg) AD0o; zeros(mg) zeros(mg,(days-1)*mg)];
% spy(AD2)
G0o = digraph(AD1o);
or_nodes = numnodes(G0o);
% initial number of edges
init_edges = numedges(G0o);
G0o = addnode(G0o,{'S*'});
G0o = addnode(G0o,{'T*'});
Ed_S0o = table([(or_nodes+1)*ones(mg,1) (1:mg)'],ones(mg,1),'VariableNames',{'EndNodes','Weight'});
Ed_T0o = table([ [(or_nodes-mg+1):or_nodes]' (or_nodes+2)*ones(mg,1)],ones(mg,1),'VariableNames',{'EndNodes','Weight'});
G0o = addedge(G0o,Ed_S0o);
G0o = addedge(G0o,Ed_T0o);
%% add edges costs and capacities 
% edge capacities are 1 except those of S* and T* with capacity equal to
% the # of MESS
% assign label nodes
G0o.Edges.Labels = (1:numedges(G0o))';
G0o.Edges.Capacities = G0o.Edges.Labels;
% figure(2)
subplot(1,2,1)
h2 =plot(G0o);
layout(h2,'layered','Direction','right','Sources', 'S*','Sinks','T*')
highlight(h2,'Edges',1:numedges(G0o),'LineWidth',1.5)
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'Visible','off')

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

% box off
labelnode(h2,[1:numnodes(G0o)],'')
% xd = get(h2, 'XData');
% yd = get(h2, 'YData');
% text(xd, yd, '', 'FontSize',20, 'FontWeight','bold')%, 'HorizontalAlignment','left', 'VerticalAlignment','middle')
% text(h2.XData(end-1)+0.1,h2.YData(end-1),'S*','fontsize',10,'HorizontalAlignment','right')
% text(h2.XData(end)+0.1,h2.YData(end),'T*','fontsize',10)

xlabel('(a)')
highlight(h2,[1:numnodes(G0o)]' ,'NodeColor','r','MarkerSize',6)
highlight(h2,'Edges',1:numedges(G0o),'LineWidth',1.5)
% highlight(h1,1:numnodes(G0o),'MarkerSize',4)
%%
print -depsc2 Fig_3_ndsplt_6.eps


