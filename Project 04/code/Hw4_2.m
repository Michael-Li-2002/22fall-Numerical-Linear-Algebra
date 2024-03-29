%% 数值代数 第4次上机作业 李佳 2100010793
% P2: 求二次多项式
% 直接点击“运行”,工作区 sol 变量处查看结果
% 输出 e1 为 QR法 的残向量2-范数, e2 为 正则化方法 的残向量2-范数
%%
t = [-1; -0.75; -0.5; 0; 0.25; 0.5; 0.75];
y = [1.00; 0.8125; 0.75; 1.00; 1.3125; 1.75; 2.3125];
A = zeros(7,3); 
for i = 1:7   % 设置系数矩阵
    A(i,1) = t(i)^2;
end
A(:,2) = t; A(:,3) = ones(7,1);
sol = LS_sol(A,y);       % 最小二乘求解
e1 = norm(A*sol - y)     % 计算残向量2-范数
[reg_sol,~] = LUcol_sol(A'*A,A'*y); % 正则化方法求解最小二乘
e2 = norm(A*reg_sol - y) % 计算残向量2-范数