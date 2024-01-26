function [U,V] = CGforStokes(U,V,P,F,G)
% �����ݶȷ����������(AU=F-BP)
% �����ֵ U,V,P, ������ F,G
% �����ȷ�� U,V
[N,~] = size(U); D = zeros(N);
k = 0; [rF,rG,~] = residue(U,V,P,F,G,D);
rho = sum(sum(rF.^2)) + sum(sum(rG.^2));
[F1,G1,~] = residue(zeros(N,N-1),zeros(N-1,N),P,F,G,D);
r0 = sqrt(sum(sum(F1.^2)) + sum(sum(G1.^2)));
while sqrt(rho)/r0 > 10^(-10)
    k = k + 1;
    if k == 1
        pF = rF; pG = rG;
    else
        beta = rho/rrho; pF = rF + beta * pF; pG = rG + beta * pG;
    end
    [wF,wG,~] = residue(pF,pG,zeros(N),zeros(N,N-1),zeros(N-1,N),D); % ��һ������
    pTw = sum(sum(pF.* wF)) + sum(sum(pG.*wG));
    alph = rho/(pTw); U = U - alph * pF; V = V - alph * pG;
    rF = rF - alph * wF; rG = rG - alph * wG;
    rrho = rho; rho = sum(sum(rF.^2)) + sum(sum(rG.^2));
end
        

