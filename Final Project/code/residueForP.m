function [D_err] = residueForP(U,V,D)
%RESIDUEFORP 仅为减少运算量 只计算P对应残量
[N,~] = size(D);
D_err = zeros(N);
if N~=2
    D_err(1,1) = D(1,1) - N*(U(1,1)+V(1,1)); 
    for j=2:N-1
        D_err(1,j) = D(1,j) - N*(U(1,j)-U(1,j-1)+V(1,j)); 
    end
    D_err(1,N) = D(1,N) - N*(-U(1,N-1)+V(1,N)); 
    for k=2:N-1
        D_err(k,1) = D(k,1) - N*(U(k,1)+V(k,1)-V(k-1,1)); 
        for j=2:N-1
            D_err(k,j) = D(k,j) - N*(U(k,j)-U(k,j-1)+V(k,j)-V(k-1,j)); 
        end
        D_err(k,N) = D(k,N) - N*(-U(k,N-1)+V(k,N)-V(k-1,N)); 
    end
    D_err(N,1) = D(N,1) - N*(U(N,1)-V(N-1,1)); 
    for j=2:N-1
        D_err(N,j) = D(N,j) - N*(U(N,j)-U(N,j-1)-V(N-1,j)); 
    end
    D_err(N,N) = D(N,N) - N*(-U(N,N-1)-V(N-1,N)); 
    
elseif N==2
    D_err(1,1) = D(1,1) - N*(U(1,1)+V(1,1)); 
    D_err(1,2) = D(1,2) - N*(-U(1,1)+V(1,2)); 
    D_err(2,1) = D(2,1) - N*(U(2,1)-V(1,1)); 
    D_err(2,2) = D(2,2) - N*(-U(2,1)-V(1,2)); 
end
