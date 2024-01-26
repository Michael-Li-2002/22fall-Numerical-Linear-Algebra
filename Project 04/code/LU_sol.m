function [b,t] = LU_sol(A,b)
%LU_SOL Guass��ȥ���ⷽ����
[~,n] = size(A);
tic;
A = LU_fac(A);
for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t = toc;