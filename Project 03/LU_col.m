function [A,u] = LU_col(A)
% LU_COL 列主元LU分解
% 输入A为进行分解的矩阵 PA=LU
% 输出A = LU 向量u
% A的上三角(包含对角线)为上三角阵U 下三角(不含对角线)为单位下三角阵L(不含对角线)
% u为行交换(u(i)<->i)
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
