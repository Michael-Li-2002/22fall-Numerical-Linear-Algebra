function [A] = LU_fac(A)
%LU_FAC Guass消去法进行LU分解
[~,n] = size(A);
for i = 1:n-1
    if A(i,i) == 0
        disp('矩阵奇异')
        break;
    else
        A(i+1:n,i) = A(i+1:n,i)/A(i,i);
        A(i+1:n, i+1:n) = A(i+1:n, i+1:n) - A(i+1:n,i) * A(i,i+1:n);
    end
end