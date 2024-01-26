function [b,t] = LUcol_sol(A,b)
%LUcol_sol 列主元Guass消去法求解 
tic;
[A,u] = LUcol_fac(A);
[~,n] = size(A);
for k = 1:n-1  % 向量进行行交换( 计算P^(-1)b )
    mid = b(k);
    b(k) = b(u(k));
    b(u(k)) = mid;
end
for i = 1:n-1  % 前代法解 Ly=b(y记录在b中)
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for j = n:-1:2 % 回代法解 Ux=y(x记录在b中)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t = toc;

