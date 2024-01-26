function [x] = LS_sol(A,b)
%LS_SOL QR�������С��������
[m,n] = size(A);
[A,d] = QR_fac(A);
for i = 1:n
    if i<m
        c = d(i)*(b(i) + A(i+1:m,i)'*b(i+1:m));
        b(i) = b(i) - c;
        b(i+1:m) = b(i+1:m) - c * A(i+1:m,i);
    end
end
for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
x = b(1:n);


