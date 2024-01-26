function [A] = LU_norm(A)
% LU_norm 直接进行LU分解
% 输入A为要进行分解的矩阵 输出A 
% A的上三角(包含对角线)为上三角阵U 下三角(不含对角线)为单位下三角阵L(不含对角线)
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

