function [U,V] = VcycleForUzawa(level,r,nu1,nu2,U,V,F,G)
% 对共轭梯度法的用Vcycle实现的预优
% 输入网格规模 2^level, Vcycle层数 r+1, 前磨光次数 nu1, 后磨光次数 nu2
% 初值 U,V, 常数项 F,G
% 输出预优后的(zk)U,V
N = 2^level; 
if r==0
    for i = 1:nu1+nu2
        [U,V] = GS(U,V,F,G); 
    end
else
    for i=1:nu1
        [U,V] = GS(U,V,F,G);
    end
    %% 残量方程
    [F_err,G_err] = residueForUV(U,V,zeros(N,N),F,G);
    
    %% 粗网格校正
    F0 = restriction(F_err,N,1); G0 = restriction(G_err,N,2); 
    [U0,V0] = VcycleForUzawa(level-1,r-1,nu1,nu2,zeros(N/2,N/2-1),zeros(N/2-1,N/2),F0,G0);
    
    %% 拉回再磨光
    U1 = prolongation2(U0,N/2,1); V1 = prolongation2(V0,N/2,2);
    U = U + U1; V = V + V1;
    for i=1:nu2
        [U,V] = GS(U,V,F,G);
    end
end
    

