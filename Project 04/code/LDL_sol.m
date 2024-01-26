function [b,t] = LDL_sol(A,b)
% LDL_sol �Ľ���ƽ���������
% �������η�����
[~,n] = size(A);
tic;
A = LDL_fac(A);
for i = 1: n-1    % ǰ������ Lz=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n      % �� Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % �ش����� L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t = toc;         % ������, ֹͣ��ʱ, ��¼�����ʱ