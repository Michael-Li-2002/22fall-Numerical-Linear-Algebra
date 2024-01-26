%% 数值代数 第2次上机作业 李佳 2100010793
% 该脚本处理 N = 256, 512
% 运行前修改"矩阵设置"中的"N = xxx", 再单击运行
% 命令行窗口输出计算用时(该脚本还输出等待时间,故建议在工作区查看结果):
% t1: 一般的 LDL^T 分解 t2: 带状Guass消去法 t3: LDL^T 带状方法
%% 矩阵设置
N = 32;     % x,y轴各进行N等分
n = (N-1)^2; % 求解的网格点数、矩阵阶数
e = ones(n,1); e1 = ones(n,1); e2 = ones(n,1);
for i = 1:N-2
    e1(n - (N-1)*i+1) = 0;
    e2((N-1)*i) = 0;
end
% 设置初始矩阵 A0
A0 = spdiags([-e, -e2, 4*e, -e1, -e],[-N+1, -1, 0, 1, N-1],n,n);
% 设置初始矩阵的只存储对角线的矩阵 B0
B0 = zeros(n,2*N-1);
B0(:,1) = -e; B0(:,N-1) = -e2; B0(:,N) = 4*e; B0(:,N+1) = -e1; B0(:,2*N-1) = -e;
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
for i = 1: n-1    % 前代法解 Ly=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % 回代法解 L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t1 = toc;         % 求解结束, 停止记时, 下一步输出计算耗时
disp(['一般的 LDL^T 方法在 N=',num2str(N),' 时，耗时 ',num2str(t1),' 秒'])

%% Method 2: 带状Guass消去法
% 带状Guass消去
B = B0; b = b0;    % 将进行操作的矩阵、向量赋初值 B0, b0
tstart = tic;      % 开始记时
for i = 1: n-(N-1) % 进行前面的带状Guass消去
    if mod(i,N-1) == 1 
        tic        % 在计算的进程中每隔 N-1 轮记一次时,用于估算剩余计算用时
    end
    B(i,1:N-1) = B(i,1:N-1)/B(i,N);
    for j = 1:N-1
        B(i+j,1+j:N-1+j) = B(i+j,1+j:N-1+j) - B(i,1:N-1)*B(i+j,N+j);
    end
    if mod(i,N-1) == 0  % N-1 轮计算结束,估算剩余用时
        disp(['Please wait for another: ',num2str(toc *(n-i)/(N-1))])
    end
end
for i = n-(N-1)+1: n-1  % 后面规模更小 进行一般的高斯消去
    B(i,i-n+N:N-1) = B(i,i-n+N:N-1)/B(i,N);
    for j = 1:n-i
        B(i+j,N+j-n+i:N-1+j) = B(i+j,N+j-n+i:N-1+j) - B(i,i-n+N:N-1)*B(i+j,N+j);
    end
end

% 解三角形方程组
for i = 1: n-N+1        % 前代法解 Ly=b 每次只需计算N-1个分量
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* B(i,N-1:-1:1)';
end
for i = n-N+2: n-1      % 正常完成后面的前代法求解
    b(i+1:n) = b(i+1:n) - b(i)* B(i,N-1:-1:N-n+i)';
end
for j = n: -1: N        % 回代法解 Ux=y 每次只需计算N-1个分量
    b(j) = b(j)/B(j,N);
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* B(j,2*N-1:-1:N+1)';
end
for j = N-1: -1: 2      % 正常完成后面的回代法求解
    b(j) = b(j)/B(j,N);
    b(1:j-1) = b(1:j-1) - b(j)* B(j,N+j-1:-1:N+1)';
end
b(1) = b(1)/B(1,N);
t2 = toc(tstart);       % 求解结束, 停止记时, 下一步输出计算耗时
disp(['带状Guass消去法在 N=',num2str(N),' 时，耗时 ',num2str(t2),' 秒'])

%% Method 3: LDL^T 带状方法
% LDL^T 带状方法分解
% 将进行操作的矩阵 B 、向量 b 赋初值 B0, b0
% 定义中间需要的向量 v, u, w
b = b0; B = B0; v = zeros(N-1,1); u = zeros(1,N-1); w = zeros(N-1,1);
tstart = tic;   % 开始记时
for j = 1:n     
    if mod(j,N-1) == 1 
        tic     % 在计算的进程中每隔 N-1 轮记一次时,用于估算剩余计算用时
    end
    % 前 N-1 列的计算
    if j <= N-1  
        for i = 1:j-1
            u(i) = B(i,N-j+i);  % 对应行向量 A(j,1:j-1)
            v(i) = B(i,N-j+i)*B(i,N);
        end
        B(j,N) = B(j,N) - u(1:j-1)* v(1:j-1); % 计算 A(j,j)
        for i = N-1:-1:1  % 计算 A(j+1:j+N-1,j)
            if i >= j 
                % 可以从 带状区域最左侧开始 累加
                for k = 1:j-1
                    B(j,i) = B(j,i) - B(k,k+i-j)*v(k);
                end
            else
                % 超出带状区域边界 从后面的某一项开始 累加
                for k = 1:i-1
                    B(j,i) = B(j,i) - B(j-i+k,k)*v(j-i+k);
                end
            end
            B(j,i) = B(j,i)/B(j,N); % 得到 A(j+N-i,j)
        end
    % 计算中间列
    elseif j <= n-(N-1)
        for i = j-N+1:j-1
            u(i-(j-N)) = B(i,N-j+i);        % 对应行向量 A(j,1:j-1)
            v(i-(j-N)) = B(i,N-j+i)*B(i,N); % 向量 v 此时只需计算(N-1)个分量
        end
        B(j,N) = B(j,N) - u(1:N-1) * v(1:N-1);  % 计算 A(j,j)
        for i = N-1:-1:1  % 计算 A(j+1:j+N-1,j)
            for k = 1:i-1
                B(j,i) = B(j,i) - B(j-i+k,k)*v(N-i+k);
            end
            B(j,i) = B(j,i)/B(j,N); % 得到 A(j+N-i,j)
        end
    % 计算后 N-1 列
    else
        for i = j-N+1:j-1
            u(i-(j-N)) = B(i,N-j+i);        % 对应行向量 A(j,1:j-1)
            v(i-(j-N)) = B(i,N-j+i)*B(i,N); % 向量 v 此时只需计算(N-1)个分量
        end
        B(j,N) = B(j,N) - u(1:N-1) * v(1:N-1);  % 计算 A(j,j)
        for i = N-1:-1:N-n+j  % 计算 A(j+1:j+n,j)
            for k = 1:i-1
                B(j,i) = B(j,i) - B(j-i+k,k)*v(N-i+k);
            end
            B(j,i) = B(j,i)/B(j,N); % 得到 A(j+N-i,j)
        end
    end
    if mod(j,N-1) == 0  % N-1 轮计算结束,估算剩余用时
        disp(['Please wait for another: ',num2str(toc *(n-j)/(N-1))])
    end
end

% 解三角形方程组
for i = 1: n-N+1      % 前代法解 Lz=b 每次只需计算N-1个分量
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* B(i,N-1:-1:1)';
end
for i = n-N+2: n-1    % 正常完成后面的前代法求解
    b(i+1:n) = b(i+1:n) - b(i)* B(i,N-1:-1:N-n+i)';
end
for i = 1: n          % 解 Dy=z
    b(i) = b(i)/ B(i,N);
end
for j = n: -1: N      % 回代法解 Dy=z, L^Tx=y 每次只需计算N-1个分量
    for i = j-N+1: j-1
        w(i-j+N) = B(i,N-j+i);
    end
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* w;
end
for j = N-1: -1: 1    % 正常完成后面的回代法求解
    for i = 1: j-1
        w(i) = B(i,N-j+i);
    end
    b(1:j-1) = b(1:j-1) - b(j)* w(1:j-1);
end
t3 = toc(tstart);     % 求解结束, 停止记时, 下一步输出计算耗时
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
