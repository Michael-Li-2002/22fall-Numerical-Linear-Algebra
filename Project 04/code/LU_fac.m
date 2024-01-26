function [A] = LU_fac(A)
%LU_FAC Guass��ȥ������LU�ֽ�
[~,n] = size(A);
for i = 1:n-1
    if A(i,i) == 0
        disp('��������')
        break;
    else
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n, i+1:n) = A(i+1:n, i+1:n) - A(i+1:n,i) * A(i,i+1:n);
    end
end