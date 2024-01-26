%% 数值代数 第3次上机作业 李佳 2100010793
% P2: 估计计算解的精度, 并与真实误差比较
% 直接点击“运行” 在工作区查看'error'变量 第一行为估计误差 第二行为真实误差
% 第三行为估计误差与真实误差的比值
% 输出图为 n=21~30时真实误差，估计误差随n变化的图线(同一个n,数值更大的是估计误差.)
%% 
error = zeros(3,26);
for n = 5:30
    A = - ones(n,n);
    A = tril(A) + 2* eye(n);
    A(:,n) = ones(n,1);
    x = randn(n,1);
    b = A*x; x0 = b;
    [A0,u] = LU_col(A);
    for k = 1:n-1  % 向量进行行交换( 计算P^(-1)b )
        mid = x0(k);
        x0(k) = x0(u(k));
        x0(u(k)) = mid;
    end
    for i = 1:n-1  % 前代法解 Ly=b(y记录在b中)
        x0(i+1:n) = x0(i+1:n) - x0(i)* A0(i+1:n,i);
    end
    for j = n:-1:2 % 回代法解 Ux=y(x记录在b中)
        x0(j) = x0(j)/A0(j,j);
        x0(1:j-1) = x0(1:j-1) - x0(j)* A0(1:j-1,j);
    end
    x0(1) = x0(1)/A0(1,1);
    
    r = b - A*x0;
    rnorm = max(abs(r)); % r 的无穷范数
    bnorm = max(abs(b)); % b 的无穷范数
    Anorm = sum(abs(A(n,:))); % A 的无穷范数(显然A的行和最大在最后一行)
    AInvnorm = InftyNorm_Inv(A); % A^-1的无穷范数
    error(1,n-4) = rnorm * Anorm * AInvnorm / bnorm; %计算解的误差估计
    error(2,n-4) = max(abs(x - x0));
    error(3,n-4) = error(1,n-4)/error(2,n-4);
end
%% 估计误差与真实误差的比较
n= 21:30;
plot(n,error(1,17:26))
hold on;
plot(n,error(2,17:26))