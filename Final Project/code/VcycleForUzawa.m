function [U,V] = VcycleForUzawa(level,r,nu1,nu2,U,V,F,G)
% �Թ����ݶȷ�����Vcycleʵ�ֵ�Ԥ��
% ���������ģ 2^level, Vcycle���� r+1, ǰĥ����� nu1, ��ĥ����� nu2
% ��ֵ U,V, ������ F,G
% ���Ԥ�ź��(zk)U,V
N = 2^level; 
if r==0
    for i = 1:nu1+nu2
        [U,V] = GS(U,V,F,G); 
    end
else
    for i=1:nu1
        [U,V] = GS(U,V,F,G);
    end
    %% ��������
    [F_err,G_err] = residueForUV(U,V,zeros(N,N),F,G);
    
    %% ������У��
    F0 = restriction(F_err,N,1); G0 = restriction(G_err,N,2); 
    [U0,V0] = VcycleForUzawa(level-1,r-1,nu1,nu2,zeros(N/2,N/2-1),zeros(N/2-1,N/2),F0,G0);
    
    %% ������ĥ��
    U1 = prolongation2(U0,N/2,1); V1 = prolongation2(V0,N/2,2);
    U = U + U1; V = V + V1;
    for i=1:nu2
        [U,V] = GS(U,V,F,G);
    end
end
    

