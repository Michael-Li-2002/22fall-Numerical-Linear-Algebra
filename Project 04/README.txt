主程序：
Hw4_1:  QR法求解ch1上机习题中的3个线性方程组并比较结果
Hw4_2:  求二次多项式最小二乘结果
Hw4_3:  求出模型中参数的最小二乘结果

函数用途：
house:  计算householder变换 输出 H=I-beta*v*v^T 中的 beta, v
QR_fac: 基于Householder变换的QR分解
QR_sol: 基于 Householder 方法的 QR 分解方法解方程组
LS_sol:   QR法求解最小二乘问题
LU_fac:  Guass消去法进行LU分解  LU_sol:  Guass消去法解方程组
LU_col:  列主元LU分解  LUcol_sol:  列主元Guass消去法求解
Cholesky_fac:  一般的Cholesky分解  Cholesky_sol:  平方根法求解
LDL_fac:   LDL^T分解  LDL_sol:  改进的平方根法求解

复现上机报告结果的方法：
Hw4_1:  直接点击“运行”
              工作区 error 变量处查看解的无穷范数误差；t 变量处查看计算时间
              第i行对应第i个方程组 第1-5列分别对应 QR法、Guass消去法、列主元Guass消去法、平方根法、改进的平方根法
Hw4_2:  直接点击“运行”,工作区 sol 变量处查看结果
              输出 e1 为 QR法 的残向量2-范数, e2 为 正则化方法 的残向量2-范数
Hw4_3:  直接点击“运行”,工作区 sol 变量处查看结果
              输出 e3 为 QR法 的残向量2-范数, e4 为 正则化方法 的残向量2-范数