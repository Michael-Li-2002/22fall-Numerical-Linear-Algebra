function [v,beta] = house(x)
% 计算 Householder 变换 输出 H=I-beta*v*v^T 中的 beta, v
% try: 规格化存储放在 QR_fac 函数中
n = length(x);
x = x/max(abs(x));      % 归一化计算防止上溢下溢
sigma = x(2:n)'*x(2:n); 
v = x;
if sigma == 0    % 若向量除第一分量外均为 0 则无需进行 Householder 变换
    beta = 0;    % 令 beta 为 0 不进行 Householder 变换
else
    alpha = sqrt(x(1)^2 + sigma); % 向量 x 的 模长
    if x(1) <= 0 % 若第一分量不大于 0 则进行一般的计算
        v(1) = x(1) - alpha;
    else         % 若第一分量大于 0 则调整计算方式
        v(1) = -sigma/(x(1)+alpha);
    end          
    beta = 2 /(sigma + v(1)^2); % 计算系数 beta
end