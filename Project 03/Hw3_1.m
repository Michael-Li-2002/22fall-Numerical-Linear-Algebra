%% 数值代数 第3次上机作业 李佳 2100010793
% P1: 估计 5-20 阶 Hilbert 矩阵的条件数
% 直接点击“运行” 在工作区查看变量'cond_num'中计算的条件数估计值
%% 设置矩阵
cond_num = zeros(1,16);
%% 估计条件数
for n = 5:20
    A = hilb(n);
    InfNormA = sum(A(1,:));
    InfNormAInv = InftyNorm_Inv(A);
    cond_num(n-4) = InfNormA * InfNormAInv;
end