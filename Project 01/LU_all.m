function [A,u,v] = LU_all(A)
% LU_all ȫ��Ԫ���Ƿֽ�
% ����AΪ���зֽ�ľ��� PAQ=LU
% ���A = LU ����u,v
% A��������(�����Խ���)Ϊ��������U ������(�����Խ���)Ϊ��λ��������L(�����Խ���)
% uΪ�н���(u(i)<->i)��vΪ�н���(v(i)<->i)
[~,n] = size(A);
u = zeros(n,1); v = zeros(n,1);
for i = 1:n-1
    mainEle = 0; p=i; q=i;
    for j = i:n
        for k = i:n
            if abs(A(j,k)) > mainEle
                p = j; q = k;
                mainEle = abs(A(j,k));
            end
        end
    end
    temp1 = A(i,1:n); A(i,1:n) = A(p,1:n); A(p,1:n) = temp1;
    temp2 = A(1:n,i); A(1:n,i) = A(1:n,q); A(1:n,q) = temp2;
    u(i) = p; v(i) = q;
    if A(i,i) ~= 0
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
    else
        disp('��������')
        break;
    end
end

