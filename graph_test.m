clc;
m = 4;
d = 7;
A = ones(m,d);
% A(1,:) = 1:d
% A(2,:) = (d+1):2*d
% A(3,:) = (2*d+1):3*d
% A(4,:) = (3*d+1):4*d
for i_m = 1:m
    A(i_m,:) = ((i_m-1)*d+1):i_m*d;
end
A
B = zeros(m,d);
for i_d = 1:d
    B(:,i_d) = ((i_d-1)*m+1):i_d*m;
end
B
disp(['Indices are',num2str(B(5:8))])
Ad = zeros(size(B))
% for i_ad = 1:m
%     for i_am = 1:d
%         
%     end
% end
%%
clc
for i_m2 = 1:m
%     for i_d2 = 1:d
        [B(i_m2)*ones(m,1) B(:,i_d2)]
%     end
end
%%
clc
for i_c = 1:d-1
    [B(:,i_c) B(:,i_c+1)]
end
%%
clc; clear all;
m = 2; d = 3
Ad_g = zeros(m*d)
for i_m = 1:m
    Ad_g(i_m,:) = [zeros(1,m) ones(1,m) zeros(1,m*d-2*m)]
end
Ad_g(3,6) = 1
Ad_g(3,5) = 1
Ad_g(4,5) = 1
Ad_g(4,6) = 1

%%
% block matrix with ones(m)
clc; clear all; close all;
m=5; d = 7;
% Ad_g =  blkdiag(ones(m), ones(m), ones(m), ones(m));
Ad_g = kron(eye(d-1),ones(m));
Ad_g = [zeros(m*d-m,m) Ad_g];
Ad_g = [Ad_g; zeros(m,m*d)]
%
G = digraph(Ad_g);
plot(G)