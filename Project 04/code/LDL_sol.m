function [b,t] = LDL_sol(A,b)
% LDL_sol 改进的平方根法求解
% 解三角形方程组
[~,n] = size(A);
tic;
A = LDL_fac(A);
for i = 1: n-1    % 前代法解 Lz=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n      % 解 Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % 回代法解 L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t = toc;         % 求解结束, 停止记时, 记录计算耗时