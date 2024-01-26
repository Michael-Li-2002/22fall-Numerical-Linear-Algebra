function [X1] = prolongation2(X,N,n)
% 提升算子2(实际使用)
% 输入向量 X, 规模 N, 类型 n (U:n=1, V:n=2, P:n=3)
% 输出提升后的向量X1
N1 = 2*N;
if n==1
    X1 = zeros(N1,N1-1);
    for i = 2:2:N1
        X1(i,1) = 0.5 * X(i/2,1); X1(i,N1-1) = 0.5 * X(i/2,N-1);
        X1(i-1,1) = 0.5 * X(i/2,1); X1(i,N1-1) = 0.5 * X(i/2,N-1);
    end
    for j=2:2:N1-2
        X1(1,j) = X(1,j/2);
        X1(N1,j) = X(N,j/2);
        for i = 2:2:N1-1
            X1(i,j) = X(i/2,j/2) ;
            X1(i+1,j) = X(i/2+1,j/2);
        end
    end
    for j=3:2:N1-3
        for i = 2:2:N1
            X1(i,j) = 0.5 * ( X(i/2,(j-1)/2) + X(i/2,(j+1)/2) );
            X1(i-1,j) = 0.5 * ( X(i/2,(j-1)/2) + X(i/2,(j+1)/2) );
        end
    end
    
elseif n==2
    X1 = prolongation2(X',N,1)';
elseif n==3
    X1 = zeros(N1);
    for i=2:2:N1
        for j=2:2:N1
            X1(i,j) = X(i/2,j/2); X1(i-1,j) = X(i/2,j/2); X1(i,j-1) = X(i/2,j/2); X1(i-1,j-1) = X(i/2,j/2);
        end
    end
end


