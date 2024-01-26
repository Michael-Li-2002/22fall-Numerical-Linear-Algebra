function [norm1] = InftyNorm_Inv(A)
% 输入矩阵A 输出矩阵A^-1的无穷范数
OK = 1; % 用于判定是否终止循环
[~,n] = size(A);
[A0,u] = LU_col(A);% 提前作列主元三角分解
x = ones(n,1)/n; B0 = A0'; % 定义初始向量 对矩阵作转置
while OK == 1 
    w = x;  
    % 求解 A^Tw=x
    for i = 1:n-1  % 前代法解 U^Ty=b(y记录在b中)
        w(i) = w(i)/B0(i,i);
        w(i+1:n) = w(i+1:n) - w(i)* B0(i+1:n,i);
    end
    w(n) = w(n)/B0(n,n);
    for j = n:-1:2 % 回代法解 L^Tx=y(x记录在b中)
        w(1:j-1) = w(1:j-1) - w(j)* B0(1:j-1,j);
    end
    for k = n-1:-1:1  % 向量进行行交换( 计算P^Tb )
        mid = w(k);
        w(k) = w(u(k));
        w(u(k)) = mid;
    end
    z = sign(w);
    % 求解 Az=v
    for k = 1:n-1  % 向量进行行交换( 计算P^(-1)b )
        mid = z(k);
        z(k) = z(u(k));
        z(u(k)) = mid;
    end
    for i = 1:n-1  % 前代法解 Ly=z(y记录在z中)
        z(i+1:n) = z(i+1:n) - z(i)* A0(i+1:n,i);
    end
    for j = n:-1:2 % 回代法解 Ux=y(x记录在z中)
        z(j) = z(j)/A0(j,j);
        z(1:j-1) = z(1:j-1) - z(j)* A0(1:j-1,j);
    end
    z(1) = z(1)/A0(1,1);
    inftynorm = 0; % z的无穷范数; 
    pos = 1; % z某个分量的模与无穷范数相等的坐标
    for index = 1:n
        if abs(z(index)) > inftynorm
            inftynorm = abs(z(index));
            pos = index;
        end
    end
    if inftynorm > z'* x % 不是最大值 继续迭代
        x = zeros(n,1);
        x(pos) = 1;
    else % 停止迭代
        norm1 = sum(abs(w));
        OK = 0;
    end
end
    

