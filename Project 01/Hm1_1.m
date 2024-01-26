%% 数值代数 第1次上机作业 李佳 2100010793
%  实现Guass消元法解线性方程组
%% 直接点击“运行”即可看到结果（命令行窗口 或 工作区）
% 数据设置
num = [2,12,24,48,84]; % 矩阵规模
error1 = zeros(3,5);   % 2范数误差 第1行：直接Guass消去 第2行：列主元 第3行：全主元
error2 = zeros(3,5);   % 无穷范数误差

% 
for index = 1:5 
    n = num(index);
    
    A0 = zeros(n,n); b0 = zeros(n,1); % A0:原始方程中的矩阵 b0:原始方程中的向量
    realSol = ones(n,1); % 精确解向量
    % 设置原始矩阵、向量的数值
    A0(1,1) = 6; A0(2,1) = 8; A0(n-1,n) = 1; A0(n,n) = 6; 
    b0(1,1) = 7; b0(n,1) = 14;
    for i = 2:n-1 
        A0(i-1,i) = 1;
        A0(i,i) = 6;
        A0(i+1,i) = 8;
        b0(i,1) = 15;
    end

    % 1 直接Guass消去法
    A = LU_norm(A0); b = b0; % LU_norm: 直接进行Guass消去法
    for i = 1:n-1  % 前代法解 Ly=b(y记录在b中)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2 % 回代法解 Ux=y(x记录在b中)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1) = b(1)/A(1,1);
        % 误差计算
    error1(1,index) = norm(b - realSol);      % 计算解与精确解误差的2范数 
    error2(1,index) = max(abs(b - realSol));  % 计算解与精确解误差的无穷范数 
    
    % 2 列主元Guass消去法
    [A,u] = LU_col(A0); b = b0; % LU_col: 列主元Guass消去法 返回A=LU 向量u(包含行交换信息)
        
    for k = 1:n-1  % 向量进行行交换( 计算P^(-1)b )
        mid = b(k);
        b(k) = b(u(k));
        b(u(k)) = mid;
    end
    for i = 1:n-1  % 前代法解 Ly=b(y记录在b中)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2 % 回代法解 Ux=y(x记录在b中)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1) = b(1)/A(1,1);
        % 误差计算
    error1(2,index) = norm(b - realSol);     % 计算解与精确解误差的2范数
    error2(2,index) = max(abs(b - realSol)); % 计算解与精确解误差的无穷范数
    
    % 3 全主元Guass消去法
    [A,u,v] = LU_all(A0); b = b0; %LU_all: 全主元Guass消去法 返回A=LU 向量u(包含行交换信息) 向量v(包含列交换信息)

    for k = 1:n-1    % 向量进行行交换( 计算P^(-1)b )
        mid = b(k);
        b(k) = b(u(k));
        b(u(k)) = mid;
    end
    for i = 1:n-1    % 前代法解 Ly=b(y记录在b中)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2   % 回代法解 Ux=y(x记录在b中)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1,1) = b(1,1)./A(1,1);
    for k = n-1:-1:1 %向量进行行交换(计算Q^(-1)b )
        mid = b(k);
        b(k) = b(v(k));
        b(v(k)) = mid;
    end
        % 计算误差
    error1(3,index) = norm(b - realSol);     % 计算解与精确解误差的2范数
    error2(3,index) = max(abs(b - realSol)); % 计算解与精确解误差的无穷范数
end
%% 输出结果
error1
error2
