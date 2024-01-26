function [A,u,v] = LU_all(A)
% LU_all 全主元三角分解
% 输入A为进行分解的矩阵 PAQ=LU
% 输出A = LU 向量u,v
% A的上三角(包含对角线)为上三角阵U 下三角(不含对角线)为单位下三角阵L(不含对角线)
% u为行交换(u(i)<->i)，v为列交换(v(i)<->i)
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
        disp('矩阵奇异')
        break;
    end
end

