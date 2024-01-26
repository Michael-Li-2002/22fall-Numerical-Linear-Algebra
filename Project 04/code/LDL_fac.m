function [A] = LDL_fac(A)
% LDL_fac: LDL^T�ֽ� ��SPD����
% ����Գ���������A, ��� A=LDL^T
% L �洢�� A ��������(�������Խ���,L�ĶԽ��߾�Ϊ1); D �洢�� A �ĶԽ���.
[~,n] = size(A);
v = zeros(n,1);
for j = 1:n
    for i = 1:j-1
        v(i) = A(j,i)*A(i,i);
    end
    A(j,j) = A(j,j) - A(j,1:j-1)*v(1:j-1);
    A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1)*v(1:j-1))/A(j,j);
end