function [A] = Cholesky_fac(A)
% Cholesky: 一般的Cholesky分解 对SPD适用
% 输入对称正定矩阵A, 输出 A=LL^T
% L 存储在 A 的下三角
[~,n] = size(A);
for k = 1:n
    A(k,k) = sqrt(A(k,k));
    A(k+1:n,k) = A(k+1:n,k)/A(k,k);
    for j=k+1:n
        A(j:n,j) = A(j:n,j)-A(j:n,k)*A(j,k);
    end
end