function [b,t] = Cholesky_sol(A,b)
%CHOLESKY_SOL 平方根法求解
tic;
A = Cholesky_fac(A);
[~,n] = size(A);
for i = 1:n-1  % 前代法解 Ly=b(y记录在b中)
    b(i) = b(i)/A(i,i);
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
b(n) = b(n)/A(n,n);
for j = n:-1:2 % 回代法解 L^Tx=y(x记录在b中)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
b(1) = b(1)/A(1,1);
t = toc;