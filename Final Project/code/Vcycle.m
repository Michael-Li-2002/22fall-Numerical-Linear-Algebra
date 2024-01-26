function [U,V,P] = Vcycle(level,r,nu1,nu2,U,V,P,F,G,D,smoother)
% Vcycle多重网格方法
% 输入网格规模 2^level, Vcycle层数 r+1, 前磨光次数 nu1, 后磨光次数 nu2
% 初值 U,V,P, 常数项 F,G,D, 磨光子 smoother
% 输出完成一次Vcycle后的 U,V,P
N = 2^level; 
if r==0
    for i = 1:nu1+nu2
        [U,V,P] = smoother(U,V,P,F,G,D,N); 
    end
else
    for i=1:nu1
        [U,V,P] = smoother(U,V,P,F,G,D,N);
    end
    %% 残量方程
    [F_err,G_err,D_err] = residue(U,V,P,F,G,D);
    
    %% 粗网格校正
    F0 = restriction(F_err,N,1); G0 = restriction(G_err,N,2); D0 = restriction(D_err,N,3);
    [U0,V0,P0] = Vcycle(level-1,r-1,nu1,nu2,zeros(N/2,N/2-1),zeros(N/2-1,N/2),zeros(N/2),F0,G0,D0,smoother);
    
    %% 拉回再磨光
    U1 = prolongation2(U0,N/2,1); V1 = prolongation2(V0,N/2,2); P1 = prolongation2(P0,N/2,3);
    U = U + U1; V = V + V1; P = P + P1;
    for i=1:nu2
        [U,V,P] = smoother(U,V,P,F,G,D,N);
    end
end
    