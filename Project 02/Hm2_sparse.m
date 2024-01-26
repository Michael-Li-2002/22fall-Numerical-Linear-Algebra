%% ��ֵ���� ��2���ϻ���ҵ ��� 2100010793
% �ýű����� N = 256, 512
% ����ǰ�޸�"��������"�е�"N = xxx", �ٵ�������
% �����д������������ʱ(�ýű�������ȴ�ʱ��,�ʽ����ڹ������鿴���):
% t1: һ��� LDL^T �ֽ� t2: ��״Guass��ȥ�� t3: LDL^T ��״����
%% ��������
N = 32;     % x,y�������N�ȷ�
n = (N-1)^2; % ��������������������
e = ones(n,1); e1 = ones(n,1); e2 = ones(n,1);
for i = 1:N-2
    e1(n - (N-1)*i+1) = 0;
    e2((N-1)*i) = 0;
end
% ���ó�ʼ���� A0
A0 = spdiags([-e, -e2, 4*e, -e1, -e],[-N+1, -1, 0, 1, N-1],n,n);
% ���ó�ʼ�����ֻ�洢�Խ��ߵľ��� B0
B0 = zeros(n,2*N-1);
B0(:,1) = -e; B0(:,N-1) = -e2; B0(:,N) = 4*e; B0(:,N+1) = -e1; B0(:,2*N-1) = -e;
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
for i = 1: n-1    % ǰ������ Ly=b
    b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
end
for i = 1: n
    b(i) = b(i)/ A(i,i);
end
for j = n: -1: 2  % �ش����� L^Tx=y
    b(1:j-1) = b(1:j-1) - b(j)* A(j,1:j-1)';
end
t1 = toc;         % ������, ֹͣ��ʱ, ��һ����������ʱ
disp(['һ��� LDL^T ������ N=',num2str(N),' ʱ����ʱ ',num2str(t1),' ��'])

%% Method 2: ��״Guass��ȥ��
% ��״Guass��ȥ
B = B0; b = b0;    % �����в����ľ�����������ֵ B0, b0
tstart = tic;      % ��ʼ��ʱ
for i = 1: n-(N-1) % ����ǰ��Ĵ�״Guass��ȥ
    if mod(i,N-1) == 1 
        tic        % �ڼ���Ľ�����ÿ�� N-1 �ּ�һ��ʱ,���ڹ���ʣ�������ʱ
    end
    B(i,1:N-1) = B(i,1:N-1)/B(i,N);
    for j = 1:N-1
        B(i+j,1+j:N-1+j) = B(i+j,1+j:N-1+j) - B(i,1:N-1)*B(i+j,N+j);
    end
    if mod(i,N-1) == 0  % N-1 �ּ������,����ʣ����ʱ
        disp(['Please wait for another: ',num2str(toc *(n-i)/(N-1))])
    end
end
for i = n-(N-1)+1: n-1  % �����ģ��С ����һ��ĸ�˹��ȥ
    B(i,i-n+N:N-1) = B(i,i-n+N:N-1)/B(i,N);
    for j = 1:n-i
        B(i+j,N+j-n+i:N-1+j) = B(i+j,N+j-n+i:N-1+j) - B(i,i-n+N:N-1)*B(i+j,N+j);
    end
end

% �������η�����
for i = 1: n-N+1        % ǰ������ Ly=b ÿ��ֻ�����N-1������
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* B(i,N-1:-1:1)';
end
for i = n-N+2: n-1      % ������ɺ����ǰ�������
    b(i+1:n) = b(i+1:n) - b(i)* B(i,N-1:-1:N-n+i)';
end
for j = n: -1: N        % �ش����� Ux=y ÿ��ֻ�����N-1������
    b(j) = b(j)/B(j,N);
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* B(j,2*N-1:-1:N+1)';
end
for j = N-1: -1: 2      % ������ɺ���Ļش������
    b(j) = b(j)/B(j,N);
    b(1:j-1) = b(1:j-1) - b(j)* B(j,N+j-1:-1:N+1)';
end
b(1) = b(1)/B(1,N);
t2 = toc(tstart);       % ������, ֹͣ��ʱ, ��һ����������ʱ
disp(['��״Guass��ȥ���� N=',num2str(N),' ʱ����ʱ ',num2str(t2),' ��'])

%% Method 3: LDL^T ��״����
% LDL^T ��״�����ֽ�
% �����в����ľ��� B ������ b ����ֵ B0, b0
% �����м���Ҫ������ v, u, w
b = b0; B = B0; v = zeros(N-1,1); u = zeros(1,N-1); w = zeros(N-1,1);
tstart = tic;   % ��ʼ��ʱ
for j = 1:n     
    if mod(j,N-1) == 1 
        tic     % �ڼ���Ľ�����ÿ�� N-1 �ּ�һ��ʱ,���ڹ���ʣ�������ʱ
    end
    % ǰ N-1 �еļ���
    if j <= N-1  
        for i = 1:j-1
            u(i) = B(i,N-j+i);  % ��Ӧ������ A(j,1:j-1)
            v(i) = B(i,N-j+i)*B(i,N);
        end
        B(j,N) = B(j,N) - u(1:j-1)* v(1:j-1); % ���� A(j,j)
        for i = N-1:-1:1  % ���� A(j+1:j+N-1,j)
            if i >= j 
                % ���Դ� ��״��������࿪ʼ �ۼ�
                for k = 1:j-1
                    B(j,i) = B(j,i) - B(k,k+i-j)*v(k);
                end
            else
                % ������״����߽� �Ӻ����ĳһ�ʼ �ۼ�
                for k = 1:i-1
                    B(j,i) = B(j,i) - B(j-i+k,k)*v(j-i+k);
                end
            end
            B(j,i) = B(j,i)/B(j,N); % �õ� A(j+N-i,j)
        end
    % �����м���
    elseif j <= n-(N-1)
        for i = j-N+1:j-1
            u(i-(j-N)) = B(i,N-j+i);        % ��Ӧ������ A(j,1:j-1)
            v(i-(j-N)) = B(i,N-j+i)*B(i,N); % ���� v ��ʱֻ�����(N-1)������
        end
        B(j,N) = B(j,N) - u(1:N-1) * v(1:N-1);  % ���� A(j,j)
        for i = N-1:-1:1  % ���� A(j+1:j+N-1,j)
            for k = 1:i-1
                B(j,i) = B(j,i) - B(j-i+k,k)*v(N-i+k);
            end
            B(j,i) = B(j,i)/B(j,N); % �õ� A(j+N-i,j)
        end
    % ����� N-1 ��
    else
        for i = j-N+1:j-1
            u(i-(j-N)) = B(i,N-j+i);        % ��Ӧ������ A(j,1:j-1)
            v(i-(j-N)) = B(i,N-j+i)*B(i,N); % ���� v ��ʱֻ�����(N-1)������
        end
        B(j,N) = B(j,N) - u(1:N-1) * v(1:N-1);  % ���� A(j,j)
        for i = N-1:-1:N-n+j  % ���� A(j+1:j+n,j)
            for k = 1:i-1
                B(j,i) = B(j,i) - B(j-i+k,k)*v(N-i+k);
            end
            B(j,i) = B(j,i)/B(j,N); % �õ� A(j+N-i,j)
        end
    end
    if mod(j,N-1) == 0  % N-1 �ּ������,����ʣ����ʱ
        disp(['Please wait for another: ',num2str(toc *(n-j)/(N-1))])
    end
end

% �������η�����
for i = 1: n-N+1      % ǰ������ Lz=b ÿ��ֻ�����N-1������
    b(i+1:i+N-1) = b(i+1:i+N-1) - b(i)* B(i,N-1:-1:1)';
end
for i = n-N+2: n-1    % ������ɺ����ǰ�������
    b(i+1:n) = b(i+1:n) - b(i)* B(i,N-1:-1:N-n+i)';
end
for i = 1: n          % �� Dy=z
    b(i) = b(i)/ B(i,N);
end
for j = n: -1: N      % �ش����� Dy=z, L^Tx=y ÿ��ֻ�����N-1������
    for i = j-N+1: j-1
        w(i-j+N) = B(i,N-j+i);
    end
    b(j-N+1:j-1) = b(j-N+1:j-1) - b(j)* w;
end
for j = N-1: -1: 1    % ������ɺ���Ļش������
    for i = 1: j-1
        w(i) = B(i,N-j+i);
    end
    b(1:j-1) = b(1:j-1) - b(j)* w(1:j-1);
end
t3 = toc(tstart);     % ������, ֹͣ��ʱ, ��һ����������ʱ
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
