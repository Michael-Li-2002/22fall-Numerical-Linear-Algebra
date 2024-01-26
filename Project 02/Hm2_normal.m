%% ��ֵ���� ��2���ϻ���ҵ ��� 2100010793
% �ýű����� N = 32, 64, 128
% ����ǰ�޸�"��������"�е�"N = xxx", �ٵ�������
% �����д������������ʱ:
% t1: һ��� LDL^T �ֽ� t2: ��״Guass��ȥ�� t3: LDL^T ��״����
%% Set matrix
N = 32;      % x,y�������N�ȷ�
n = (N-1)^2; %��������������������
e = ones(n,1); e1 = ones(n,1); e2 = ones(n,1);
for i = 1:N-2
    e1(n - (N-1)*i+1) = 0;
    e2((N-1)*i) = 0;
end
% ���ó�ʼ���� A0
A0 = spdiags([-e, -e2, 4*e, -e1, -e],[-N+1, -1, 0, 1, N-1],n,n);
A0 = full(A0);
% �������� b0
b0 = zeros(n,1);
for i = 1: N-1 
    for j = 1: N-1
        b0((i-1)*(N-1)+j) = sin(pi*i/N)*sin(pi*j/N)/N^2;
    end
end

%% Method 1: һ��� LDL^T �ֽ�
% һ��� LDL^T �ֽ�
b = b0;           % ���� b ����ֵ b0
tic;              % ��ʼ��ʱ
A = LDL_norm(A0); % ���� LDL^T �ֽ�, ����ھ��� A ��

% �������η�����
for i = 1: n-1    % ǰ������ Lz=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n      % �� Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % �ش����� L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t1 = toc;         % ������, ֹͣ��ʱ, ��һ����������ʱ
disp(['һ��� LDL^T ������ N=',num2str(N),' ʱ����ʱ ',num2str(t1),' ��'])

%% Method 2: ��״Guass��ȥ��
% ��״Guass��ȥ
A = A0; b = b0;        % �����в����ľ�����������ֵ A0, b0
tic;                   % ��ʼ��ʱ
for i = 1: n-(N-1)     % ����ǰ��Ĵ�״Guass��ȥ
    A(i+1:i+N-1,i) = A(i+1:i+N-1,i)/A(i,i);
    A(i+1:i+N-1,i+1:i+N-1) = A(i+1:i+N-1,i+1:i+N-1) - A(i+1:i+N-1,i)*A(i,i+1:i+N-1);
end
for i = n-(N-1)+1: n-1 % �����ģ��С ����һ��ĸ�˹��ȥ
    A(i+1:n,i) = A(i+1:n,i)/A(i,i);
    A(i+1:n,i+1:n) = A(i+1:n,i+1:n) - A(i+1:n,i)*A(i,i+1:n);
end

% �������η�����
for i = 1: n-N+1       % ǰ������ Ly=b ÿ��ֻ�����N-1������
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* A(i+1:i+N-1,i);
end
for i = n-N+2: n-1     % ������ɺ����ǰ�������
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for j = n: -1: N       % �ش����� Ux=y ÿ��ֻ�����N-1������
    b(j) = b(j)/A(j,j);
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* A(j-N+1:j-1,j);
end
for j = N-1: -1: 2     % ������ɺ���Ļش������
    b(j) = b(j)/A(j,j);
    b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
end
b(1) = b(1)/A(1,1);
t2 = toc;              % ������, ֹͣ��ʱ, ��һ����������ʱ
disp(['��״Guass��ȥ���� N=',num2str(N),' ʱ����ʱ ',num2str(t2),' ��'])
%% Method 3: LDL^T ��״����
% LDL^T ��״�����ֽ�
% �����в����ľ��� A ������ b ����ֵ A0, b0, �����м���Ҫ������ v
A = A0; b = b0; v = zeros(N-1,1); 
tic              % ��ʼ��ʱ
for j = 1:n        
    % ǰ N-1 �еļ���
    if j <= N-1  
        for i = 1:j-1
            v(i) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,1:j-1)* v(1:j-1);
        % ֻ�����L�ڸ����е� N-1 ��
        A(j+1:j+N-1,j) = (A(j+1:j+N-1,j) - A(j+1:j+N-1,1:j-1)*v(1:j-1))/A(j,j);
    % ��(N-1)�еļ���
    elseif j > n-(N-1)  
        for i = j-N+1:j-1
            % ���� v ��ʱֻ�����(N-1)������
            v(i-(j-N)) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,j-N+1:j-1) * v(1:N-1);
        A(j+1:n,j) = (A(j+1:n,j) - A(j+1:n,j-N+1:j-1)*v(1:N-1))/A(j,j);
    % �м��еļ���
    else
        for i = j-N+1:j-1 
            % ���� v ��ʱֻ����� N-1 ������
            v(i-(j-N)) = A(j,i)*A(i,i);
        end
        A(j,j) = A(j,j) - A(j,j-N+1:j-1) * v(1:N-1);
        % ֻ����� L �ڸ����е� N-1 ��
        A(j+1:j+N-1,j) = (A(j+1:j+N-1,j) - A(j+1:j+N-1,j-N+1:j-1)*v(1:N-1))/A(j,j);
    end
end

% �������η�����
for i = 1: n-N+1      % ǰ������ Lz=b ÿ��ֻ�����N-1������
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* A(i+1:i+N-1,i);
end
for i = n-N+2: n-1    % ������ɺ����ǰ�������
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n          % �� Dy=z
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: N      % �ش����� Dy=z, L^Tx=y ÿ��ֻ�����N-1������
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* A(j,j-N+1:j-1)';
end
for j = N-1: -1: 2    % ������ɺ���Ļش������
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
b(1) = b(1)/A(1,1);
t3 = toc;             % ������, ֹͣ��ʱ, ��һ����������ʱ
disp(['LDL^T ��״������ N=',num2str(N),' ʱ,��ʱ ',num2str(t3),' ��'])

%% ���ӻ��������ڴ���ֱ�ۼ�����ֵ���׼ȷ��
Z = zeros(N+1,N+1);
for i = 1:N+1
    if i == 1 || i == N+1
        Z(i,j) = 0;
    else
        for j = 1:N+1
            if j == 1 || j == N+1
                Z(i,j) = 0;
            else
                Z(i,j) = b(i-1+(j-2)*(N-1));
            end
        end
    end
end
X = 0:1/N:1; Y = 0:1/N:1;
meshc(X,Y,Z);
