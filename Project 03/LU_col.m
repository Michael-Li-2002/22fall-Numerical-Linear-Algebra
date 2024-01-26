function [A,u] = LU_col(A)
% LU_COL ����ԪLU�ֽ�
% ����AΪ���зֽ�ľ��� PA=LU
% ���A = LU ����u
% A��������(�����Խ���)Ϊ��������U ������(�����Խ���)Ϊ��λ��������L(�����Խ���)
% uΪ�н���(u(i)<->i)
[~,n] = size(A);
u = zeros(n,1);
for k = 1:n-1
    p = k;
    main_value = 0;
    for j = k:n
        if abs(A(j,k)) > main_value
            main_value = abs(A(j,k));
            p = j;
        end
    end
    u(k) = p;
    temp = A(k,1:n);
    A(k,1:n) = A(p,1:n);
    A(p,1:n) = temp;
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    A(k+1:n,k+1:n) = A(k+1:n,k+1:n) - A(k+1:n,k)*A(k,k+1:n);
end
