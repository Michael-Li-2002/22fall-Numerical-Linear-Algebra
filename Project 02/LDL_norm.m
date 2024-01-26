function [A] = LDL_norm(A)
% LDL_norm: 一般的LDL^T算法 对SPD适用
% 输入对称正定矩阵A, 输出 A=LDL^T
% L 存储在 A 的下三角(不包括对角线,L的对角线均为1); D 存储在 A 的对角线.
[~,n] = size(A);
v = zeros(n,1);
for j = 1:n
    for i = 1:j-1
        v(i) = A(j,i)*A(i,i);
    end
    A(j,j) = A(j,j) - A(j,1:j-1)*v(1:j-1);
    A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,1:j-1)*v(1:j-1))/A(j,j);
end

