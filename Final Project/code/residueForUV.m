function [F_err,G_err] = residueForUV(U,V,P,F,G)
%RESIDUEFORUV 仅为减少运算量 只计算UV对应残量
[N,~] = size(F);
F_err = zeros(N,N-1); G_err = zeros(N-1,N);
if N~=2
    F_err(1,1) = F(1,1) - (3*N^2*U(1,1) - N^2*U(1,2) - N^2*U(2,1) + N* (P(1,2) - P(1,1)));
    for j=2:N-2
        F_err(1,j) = F(1,j) - (3*N^2*U(1,j) - N^2*U(1,j+1) -N^2*U(1,j-1) - N^2*U(2,j) + N* (P(1,j+1) - P(1,j)));
    end
    F_err(1,N-1) = F(1,N-1) - (3*N^2*U(1,N-1) - N^2*U(1,N-2) - N^2*U(2,N-1) + N* (P(1,N) - P(1,N-1)));
    for k=2:N-1
        F_err(k,1) = F(k,1) - (4*N^2*U(k,1) - N^2*U(k,2) - N^2*U(k+1,1) - N^2*U(k-1,1) + N* (P(k,2) - P(k,1)));
        for j=2:N-2
            F_err(k,j) = F(k,j) - (4*N^2*U(k,j) - N^2*U(k,j+1) -N^2*U(k,j-1) - N^2*U(k+1,j) - N^2*U(k-1,j) + N* (P(k,j+1) - P(k,j)));
        end
        F_err(k,N-1) = F(k,N-1) - (4*N^2*U(k,N-1) - N^2*U(k,N-2) - N^2*U(k+1,N-1) - N^2*U(k-1,N-1) + N* (P(k,N) - P(k,N-1)));
    end
    F_err(N,1) = F(N,1) - (3*N^2*U(N,1) - N^2*U(N,2) - N^2*U(N-1,1) + N* (P(N,2) - P(N,1)));
    for j=2:N-2
        F_err(N,j) = F(N,j) - (3*N^2*U(N,j) - N^2*U(N,j+1) -N^2*U(N,j-1) - N^2*U(N-1,j) + N* (P(N,j+1) - P(N,j)));
    end
    F_err(N,N-1) = F(N,N-1) - (3*N^2*U(N,N-1) - N^2*U(N,N-2) - N^2*U(N-1,N-1) + N* (P(N,N) - P(N,N-1)));

    G_err(1,1) = G(1,1) - (3*N^2*V(1,1) - N^2*V(1,2) - N^2*V(2,1) + N* (P(2,1) - P(1,1)));
    for j=2:N-1
        G_err(1,j) = G(1,j) - (4*N^2* V(1,j) - N^2*V(1,j+1) -N^2*V(1,j-1) - N^2*V(2,j) + N* (P(2,j) - P(1,j)));
    end
    G_err(1,N) = G(1,N) - (3*N^2*V(1,N) - N^2*V(1,N-1) - N^2*V(2,N) + N* (P(2,N) - P(1,N)));
    for k=2:N-2
        G_err(k,1) = G(k,1) - (3*N^2*V(k,1) - N^2*V(k,2) - N^2*V(k+1,1) - N^2*V(k-1,1) + N* (P(k+1,1) - P(k,1)));
        for j=2:N-1
            G_err(k,j) = G(k,j) - (4*N^2*V(k,j) - N^2*V(k,j+1) -N^2*V(k,j-1) - N^2*V(k+1,j) - N^2*V(k-1,j) + N* (P(k+1,j) - P(k,j)));
        end
        G_err(k,N) = G(k,N) - (3*N^2*V(k,N) - N^2*V(k,N-1) - N^2*V(k+1,N) - N^2*V(k-1,N) + N* (P(k+1,N) - P(k,N)));
    end
    G_err(N-1,1) = G(N-1,1) - (3*N^2*V(N-1,1) - N^2*V(N-1,2) - N^2*V(N-2,1) + N* (P(N,1) - P(N-1,1)));
    for j=2:N-1
        G_err(N-1,j) = G(N-1,j) - (4*N^2*V(N-1,j) - N^2*V(N-1,j+1) -N^2*V(N-1,j-1) - N^2*V(N-2,j) + N* (P(N,j) - P(N-1,j)));
    end
    G_err(N-1,N) = G(N-1,N) - (3*N^2*V(N-1,N) - N^2*V(N-1,N-1) - N^2*V(N-2,N) + N* (P(N,N) - P(N-1,N))); 
    
elseif N==2
    F_err(1,1) = F(1,1) - (3*N^2*U(1,1) - N^2*U(2,1) + N* (P(1,2) - P(1,1)));
    F_err(2,1) = F(2,1) - (3*N^2*U(2,1) - N^2*U(1,1) + N* (P(2,2) - P(2,1)));
    G_err(1,1) = G(1,1) - (3*N^2*V(1,1) - N^2*V(1,2) + N* (P(2,1) - P(1,1)));
    G_err(1,2) = G(1,2) - (3*N^2*V(1,2) - N^2*V(1,1) + N* (P(2,2) - P(1,2)));
end

