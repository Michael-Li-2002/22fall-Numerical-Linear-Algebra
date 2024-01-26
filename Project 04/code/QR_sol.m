function [b,t] = QR_sol(A,b)
% ���� Householder ������ QR �ֽⷽ���ⷽ����
% ���������Ľ�
tic;
[A,d] = QR_fac(A);   % �Ծ������ QR �ֽ�
[~,n] = size(A);
for i = 1:n-1        % �������� Q^Tb
    c = d(i)*(b(i) + A(i+1:n,i)'*b(i+1:n)); % ���� beta(v^Tb)
    b(i) = b(i) - c; % �ֱ���� b(i) �� b(i+1:n)
    b(i+1:n) = b(i+1:n) - c * A(i+1:n,i);
end
for j = n:-1:2       % �ش����� Ux=y
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t = toc;

