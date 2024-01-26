%% 数值代数 第4次上机作业 李佳 2100010793
% P1: 求解ch1上机习题中的3个线性方程组并比较结果
% 直接点击“运行”
% 工作区 error 变量处查看解的无穷范数误差；t 变量处查看计算时间
% 第i行对应第i个方程组 第1-5列分别对应 QR法、Guass消去法、列主元Guass消去法、平方根法、改进的平方根法
%% 矩阵设置
% 方程组(1)
e = ones(84,1); b1 = 15* e;
A1 = full(spdiags([8*e,6*e,e],-1:1,84,84));
b1(1) = 7; b1(84) = 14;
realSol1 = ones(84,1); % 精确解向量
% 方程组(2)
e = ones(100,1); b2 = 12* e;
A2 = full(spdiags([e,10*e,e], -1:1, 100, 100));
b2(1) = 11; b2(100) = 11; 
realSol2 = ones(100,1); % 精确解向量
% 方程组(3)
A3 = hilb(40); b3 = zeros(40,1);
for i=1:40
    b3(i) = sum(A3(i,:));
end
realSol3 = ones(40,1); % 精确解向量
%
error = zeros(3,5); t = zeros(3,5);
%% QR分解法
[x1,t1] = QR_sol(A1,b1); error(1,1)=norm(x1-realSol1,inf);t(1,1) = t1;
[x2,t2] = QR_sol(A2,b2); error(2,1)=norm(x2-realSol2,inf);t(2,1) = t2;
[x3,t3] = QR_sol(A3,b3); error(3,1)=norm(x3-realSol3,inf);t(3,1) = t3;
%% Guass消去法
[x1,t1] = LU_sol(A1,b1); error(1,2)=norm(x1-realSol1,inf);t(1,1) = t1;
%% 列主元Guass消去法
[x1,t1] = LUcol_sol(A1,b1); error(1,3)=norm(x1-realSol1,inf);t(1,3) = t1;
%% Cholesky法
[x2,t2] = Cholesky_sol(A2,b2); error(2,4)=norm(x2-realSol2,inf);t(2,4) = t2;
[x3,t3] = Cholesky_sol(A3,b3); error(3,4)=norm(x3-realSol3,inf);t(3,4) = t3;
%% LDL^T法
[x2,t2] = LDL_sol(A2,b2); error(2,5)=norm(x2-realSol2,inf);t(2,5) = t2;
[x3,t3] = LDL_sol(A3,b3); error(3,5)=norm(x3-realSol3,inf);t(3,5) = t3;