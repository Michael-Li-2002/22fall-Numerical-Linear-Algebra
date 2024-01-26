function [b,t] = Cholesky_sol(A,b)
%CHOLESKY_SOL ƽ���������
tic;
A = Cholesky_fac(A);
[~,n] = size(A);
for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
    b(i) = b(i)/A(i,i);
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
b(n) = b(n)/A(n,n);
for j = n:-1:2 % �ش����� L^Tx=y(x��¼��b��)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
b(1) = b(1)/A(1,1);
t = toc;