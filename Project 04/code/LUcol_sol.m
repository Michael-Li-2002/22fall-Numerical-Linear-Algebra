function [b,t] = LUcol_sol(A,b)
%LUcol_sol ����ԪGuass��ȥ����� 
tic;
[A,u] = LUcol_fac(A);
[~,n] = size(A);
for k = 1:n-1  % ���������н���( ����P^(-1)b )
    mid = b(k);
    b(k) = b(u(k));
    b(u(k)) = mid;
end
for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t = toc;

