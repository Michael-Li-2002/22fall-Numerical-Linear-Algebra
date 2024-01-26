function [b,t] = QR_sol(A,b)
% 基于 Householder 方法的 QR 分解方法解方程组
% 输出方程组的解
tic;
[A,d] = QR_fac(A);   % 对矩阵进行 QR 分解
[~,n] = size(A);
for i = 1:n-1        % 计算向量 Q^Tb
    c = d(i)*(b(i) + A(i+1:n,i)'*b(i+1:n)); % 计算 beta(v^Tb)
    b(i) = b(i) - c; % 分别计算 b(i) 和 b(i+1:n)
    b(i+1:n) = b(i+1:n) - c * A(i+1:n,i);
end
for j = n:-1:2       % 回代法解 Ux=y
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t = toc;

