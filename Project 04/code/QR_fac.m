function [A, d] = QR_fac(A)
% 基于Householder变换的QR分解
% try: 规格化存储放在这里
[m,n] = size(A);
d = zeros(n,1); % 存储 Householder 变换中对应的系数
for j = 1:n
    if j < m    
        [v,beta] = house(A(j:m,j)); % 计算应进行的 Householder 变换
        % 对子矩阵进行 Householder 变换
        % 将矩阵矩阵运算化为矩阵向量运算与向量向量运算
        A(j:m,j:n) = A(j:m, j:n) - (beta * v) * (v' * A(j:m, j:n)); 
        d(j) = beta * v(1)^2;         % 存储规格化后的系数 beta
        A(j+1:m,j) = v(2:m-j+1)/v(1); % 存储规格化后的向量后面的分量
    end
end

