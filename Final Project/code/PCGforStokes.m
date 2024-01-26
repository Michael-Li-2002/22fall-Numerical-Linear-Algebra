function [U,V] = PCGforStokes(U,V,P,F,G,tau,level,r,nu1,nu2)
% Ԥ�Ź����ݶȷ���������� ��Vcycle��Ԥ��
% �����ֵ U,V,P, ������ F,G, ���ƽ⾫��Ҫ�� tau, �����ģ 2^level
% Vcycle���� r+1, ǰĥ����� nu1, ��ĥ����� nu2
% ����ﵽ���ƽ⾫��Ҫ��Ľ� U,V
[N,~] = size(U); D = zeros(N);
k = 0; [rF,rG] = residueForUV(U,V,P,F,G);
res = sqrt(sum(sum(rF.^2)) + sum(sum(rG.^2)));
BTU = residueForP(U,V,D); r0 = sqrt(sum(sum(BTU.^2)));
if r0 == 0
    r0 = 1;
end
while res > tau*r0
    [zF,zG] = VcycleForUzawa(level,r,nu1,nu2,zeros(N,N-1),zeros(N-1,N),rF,rG);
    k = k + 1;
    if k == 1
        pF = zF; pG = zG; rho = sum(sum(rF.*zF)) + sum(sum(rG.*zG));
    else
        rrho = rho; rho = sum(sum(rF.*zF)) + sum(sum(rG.*zG));
        beta = rho/rrho; pF = zF + beta * pF; pG = zG + beta * pG;
    end
    [wF,wG] = residueForUV(pF,pG,zeros(N),zeros(N,N-1),zeros(N-1,N)); % ��һ������
    pTw = sum(sum(pF.* wF)) + sum(sum(pG.* wG));
    alph = rho/(pTw); U = U - alph * pF; V = V - alph * pG;
    rF = rF - alph * wF; rG = rG - alph * wG;
    res = sqrt(sum(sum(rF.^2)) + sum(sum(rG.^2)));
    %disp(['PCG complete ',num2str(k),' times, with residue ',num2str(res)])
end

