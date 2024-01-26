%% ��ֵ���� ��3���ϻ���ҵ ��� 2100010793
% P2: ���Ƽ����ľ���, ������ʵ���Ƚ�
% ֱ�ӵ�������С� �ڹ������鿴'error'���� ��һ��Ϊ������� �ڶ���Ϊ��ʵ���
% ������Ϊ�����������ʵ���ı�ֵ
% ���ͼΪ n=21~30ʱ��ʵ�����������n�仯��ͼ��(ͬһ��n,��ֵ������ǹ������.)
%% 
error = zeros(3,26);
for n = 5:30
    A = - ones(n,n);
    A = tril(A) + 2* eye(n);
    A(:,n) = ones(n,1);
    x = randn(n,1);
    b = A*x; x0 = b;
    [A0,u] = LU_col(A);
    for k = 1:n-1  % ���������н���( ����P^(-1)b )
        mid = x0(k);
        x0(k) = x0(u(k));
        x0(u(k)) = mid;
    end
    for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
        x0(i+1:n) = x0(i+1:n) - x0(i)* A0(i+1:n,i);
    end
    for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
        x0(j) = x0(j)/A0(j,j);
        x0(1:j-1) = x0(1:j-1) - x0(j)* A0(1:j-1,j);
    end
    x0(1) = x0(1)/A0(1,1);
    
    r = b - A*x0;
    rnorm = max(abs(r)); % r �������
    bnorm = max(abs(b)); % b �������
    Anorm = sum(abs(A(n,:))); % A �������(��ȻA���к���������һ��)
    AInvnorm = InftyNorm_Inv(A); % A^-1�������
    error(1,n-4) = rnorm * Anorm * AInvnorm / bnorm; %������������
    error(2,n-4) = max(abs(x - x0));
    error(3,n-4) = error(1,n-4)/error(2,n-4);
end
%% �����������ʵ���ıȽ�
n= 21:30;
plot(n,error(1,17:26))
hold on;
plot(n,error(2,17:26))