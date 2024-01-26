主程序：
StokesDGS:  第一问(以DGS为磨光子的多重网格方法求解Stokes方程）
StokesExactUzawa:  第二问(Exact Uzawa方法求解Stokes方程）
StokesInexactUzawa:  第三问(Inexact Uzawa方法求解Stokes方程）

函数用途：
f,g: 算例中的函数f(x,y),g(x,y)
residue: 求方程残量
residueForP: 只求散度项残量
residueForUV: 只求速度项残量
DGS: 第一问磨光子
prolongation: PPT上V-cycle的提升算子
prolongation2: V-cycle实际使用的提升算子
restriction: PPT上V-cycle的限制算子
Vcycle: 第一问的V-cycle多重网格方法
CGforStokes: 共轭梯度法求解子问题AU=F-BP
GS: 第三问中预条件子V-cycle的磨光子
PCGforStokes: 第三问预优共轭梯度法
VcycleForUzawa: 第三问的预条件子Vcycle

复现上机报告结果的方法：
各主程序中, 在设置参数处调整合适的参数后单击运行即可. 