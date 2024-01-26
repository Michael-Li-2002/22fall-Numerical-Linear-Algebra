%% 数值代数大作业 第二问(Exact Uzawa方法求解Stokes方程）
% 设置参数处调整合适的参数后单击运行即可
%% 设置参数
level = 7; alpha = 1; % 网格规模 N=2^level; 迭代法参数 alpha
%% 设置方程组常数项与真解
N = 2^level; h = 1/N; % 网格单位长度 h
F = zeros(N,N-1); G = zeros(N-1,N); D = zeros(N);
for j=1:N-1 % 定义线性方程组的常数项
    F(1,j) = f(j*h,0.5*h) - 2*pi*(1-cos(2*pi*j*h))*N;
    F(N,j) = f(j*h,1-0.5*h) + 2*pi*(1-cos(2*pi*j*h))*N;
    G(j,1) = g(0.5*h,j*h) + 2*pi*(1-cos(2*pi*j*h))*N;
    G(j,N) = g(1-0.5*h,j*h) - 2*pi*(1-cos(2*pi*j*h))*N;
    for i=2:N-1
        F(i,j) = f(j*h,(i-0.5)*h);
        G(j,i) = g((i-0.5)*h,j*h);
    end
end
real_U = zeros(N,N-1); real_V = zeros(N-1,N); real_P = zeros(N); 
for i=1:N % 定义真解
    for j = 1: N-1
        real_U(i,j) = (1-cos(2*pi*j*h))*sin(2*pi*(i-0.5)*h);
        real_V(j,i) = -(1-cos(2*pi*j*h))*sin(2*pi*(i-0.5)*h);
    end
    for j = 1:N
        real_P(i,j) = (j-0.5)^3*h^3/3-1/12;
    end
end
%% 求解Stokes方程
U = zeros(N,N-1); V = zeros(N-1,N); P = zeros(N);
[F_err,G_err,D_err] = residue(U,V,P,F,G,D); % 求初始残量
r0 = sqrt(sum(sum(F_err.^2))+sum(sum(G_err.^2))+sum(sum(D_err.^2))); % 初始残量2范数 r0
res = r0; cycle = 0; % 当前残量2范数 res; 外循环次数 cycle
while res/r0 > 10^(-8)
    cycle = cycle + 1;
    [U,V] = CGforStokes(U,V,P,F,G); % 共轭梯度法求解子问题
    [~,~,D_err] = residue(U,V,P,F,G,D); % 计算B^TU
    P = P + alpha * D_err; % 更新压力
    [F_err,G_err,D_err] = residue(U,V,P,F,G,D); % 计算当前残量
    res = sqrt(sum(sum(F_err.^2))+sum(sum(G_err.^2))+sum(sum(D_err.^2)));
    disp(['ExactUzawa complete ',num2str(cycle),' times, with residue ',num2str(res)])
end
%% 误差计算
error = h*sqrt(sum(sum((real_U-U).^2)) + sum(sum((real_V-V).^2)));
disp(['Numerical error: ',num2str(error),'.'])