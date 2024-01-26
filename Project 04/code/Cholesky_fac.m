function [A] = Cholesky_fac(A)
% Cholesky: һ���Cholesky�ֽ� ��SPD����
% ����Գ���������A, ��� A=LL^T
% L �洢�� A ��������
[~,n] = size(A);
for k = 1:n
    A(k,k) = sqrt(A(k,k));
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    for j=k+1:n
        A(j:n,j) = A(j:n,j)-A(j:n,k)*A(j,k);
    end
end