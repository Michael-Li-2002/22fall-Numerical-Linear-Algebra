function [X0] = restriction(X,N,n)
% 限制算子(课件上)
% 输入向量 X, 规模 N, 类型 n (U:n=1, V:n=2, P:n=3)
% 输出限制后的向量X0
N0 = N/2;
if n==1 % F
    X0 = zeros(N0,N0-1);
    for i=1:N0
        for j=1:N0-1
            X0(i,j) = 0.125*(X(2*i-1,2*j-1) + X(2*i-1,2*j+1) + X(2*i,2*j-1) + X(2*i,2*j+1)) + 0.25*(X(2*i-1,2*j) + X(2*i,2*j));
        end
    end
elseif n==2 % G
    X0 = restriction(X',N,1)';
elseif n==3 % D
    X0 = zeros(N0);
    for i=1:N0
        for j=1:N0
            X0(i,j) = 0.25*(X(2*i-1,2*j-1) + X(2*i-1,2*j) + X(2*i,2*j-1) + X(2*i,2*j));
        end
    end
end

