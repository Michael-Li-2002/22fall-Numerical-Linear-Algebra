%% 数值代数 第2次上机作业 李佳 2100010793
% 该脚本处理 N = 32, 64, 128
% 运行前修改"矩阵设置"中的"N = xxx", 再单击运行
% 命令行窗口输出计算用时:
% t1: 一般的 LDL^T 分解 t2: 带状Guass消去法 t3: LDL^T 带状方法
%% Set matrix
N = 32;      % x,y轴各进行N等分
n = (N-1)^2; %求解的网格点数、矩阵阶数
e = ones(n,1); e1 = ones(n,1); e2 = ones(n,1);
for i = 1:N-2
    e1(n - (N-1)*i+1) = 0;
    e2((N-1)*i) = 0;
end
% 设置初始矩阵 A0
A0 = spdiags([-e, -e2, 4*e, -e1, -e],[-N+1, -1, 0, 1, N-1],n,n);
A0 = full(A0);
% 设置向量 b0
b0 = zeros(n,1);
for i = 1: N-1 
    for j = 1: N-1
        b0((i-1)*(N-1)+j) = sin(pi*i/N)*sin(pi*j/N)/N^2;
    end
end

%% Method 1: 一般的 LDL^T 分解
% 一般的 LDL^T 分解
b = b0;           % 向量 b 赋初值 b0
tic;              % 开始记时
A = LDL_norm(A0); % 进行 LDL^T 分解, 输出在矩阵 A 中

% 解三角形方程组
for i = 1: n-1    % 前代法解 Lz=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n      % 解 Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % 回代法解 L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t1 = toc;         % 求解结束, 停止记时, 下一步输出计算耗时
disp(['一般的 LDL^T 方法在 N=',num2str(N),' 时，耗时 ',num2str(t1),' 秒'])

%% Method 2: 带状Guass消去法
% 带状Guass消去
A = A0; b = b0;        % 将进行操作的矩阵、向量赋初值 A0, b0
tic;                   % 开始记时
for i = 1: n-(N-1)     % 进行前面的带状Guass消去
    A(i+1:i+N-1,i) = A(i+1:i+N-1,i)/A(i,i);
    A(i+1:i+N-1,i+1:i+N-1) = A(i+1:i+N-1,i+1:i+N-1) - A(i+1:i+N-1,i)*A(i,i+1:i+N-1);
end
for i = n-(N-1)+1: n-1 % 后面规模更小 进行一般的高斯消去
    A(i+1:n,i) = A(i+1:n,i)/A(i,i);
    A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
end

% 解三角形方程组
for i = 1: n-N+1       % 前代法解 Ly=b 每次只需计算N-1个分量
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* A(i+1:i+N-1,i);
end
for i = n-N+2: n-1     % 正常完成后面的前代法求解
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for j = n: -1: N       % 回代法解 Ux=y 每次只需计算N-1个分量
    b(j) = b(j)/A(j,j);
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* A(j-N+1:j-1,j);
end
for j = N-1: -1: 2     % 正常完成后面的回代法求解
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t2 = toc;              % 求解结束, 停止记时, 下一步输出计算耗时
disp(['带状Guass消去法在 N=',num2str(N),' 时，耗时 ',num2str(t2),' 秒'])
%% Method 3: LDL^T 带状方法
% LDL^T 带状方法分解
% 将进行操作的矩阵 A 、向量 b 赋初值 A0, b0, 定义中间需要的向量 v
A = A0; b = b0; v = zeros(N-1,1); 
tic              % 开始计时
for j = 1:n        
    % 前 N-1 列的计算
    if j <= N-1  
        for i = 1:j-1
            v(i) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,1:j-1)* v(1:j-1);
        % 只需计算L在该列中的 N-1 行
        A(j+1:j+N-1,j) = (A(j+1:j+N-1,j) - A(j+1:j+N-1,1:j-1)*v(1:j-1))/A(j,j);
    % 后(N-1)列的计算
    elseif j > n-(N-1)  
        for i = j-N+1:j-1
            % 向量 v 此时只需计算(N-1)个分量
            v(i-(j-N)) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,j-N+1:j-1) * v(1:N-1);
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,j-N+1:j-1)*v(1:N-1))/A(j,j);
    % 中间列的计算
    else
        for i = j-N+1:j-1 
            % 向量 v 此时只需计算 N-1 个分量
            v(i-(j-N)) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,j-N+1:j-1) * v(1:N-1);
        % 只需计算 L 在该列中的 N-1 行
        A(j+1:j+N-1,j) = (A(j+1:j+N-1,j) - A(j+1:j+N-1,j-N+1:j-1)*v(1:N-1))/A(j,j);
    end
end

% 解三角形方程组
for i = 1: n-N+1      % 前代法解 Lz=b 每次只需计算N-1个分量
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* A(i+1:i+N-1,i);
end
for i = n-N+2: n-1    % 正常完成后面的前代法求解
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n          % 解 Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: N      % 回代法解 Dy=z, L^Tx=y 每次只需计算N-1个分量
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* A(j,j-N+1:j-1)';
end
for j = N-1: -1: 2    % 正常完成后面的回代法求解
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
b(1) = b(1)/A(1,1);
t3 = toc;             % 求解结束, 停止记时, 下一步输出计算耗时
disp(['LDL^T 带状方法在 N=',num2str(N),' 时,耗时 ',num2str(t3),' 秒'])

%% 可视化：可用于粗略直观检验数值解的准确性
Z = zeros(N+1,N+1);
for i = 1:N+1
    if i == 1 || i == N+1
        Z(i,j) = 0;
    else
        for j = 1:N+1
            if j == 1 || j == N+1
                Z(i,j) = 0;
            else
                Z(i,j) = b(i-1+(j-2)*(N-1));
            end
        end
    end
end
X = 0:1/N:1; Y = 0:1/N:1;
meshc(X,Y,Z);
