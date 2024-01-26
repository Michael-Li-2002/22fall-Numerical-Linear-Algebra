%% ��ֵ���� ��1���ϻ���ҵ ��� 2100010793
%  ʵ��Guass��Ԫ�������Է�����
%% ֱ�ӵ�������С����ɿ�������������д��� �� ��������
% ��������
num = [2,12,24,48,84]; % �����ģ
error1 = zeros(3,5);   % 2������� ��1�У�ֱ��Guass��ȥ ��2�У�����Ԫ ��3�У�ȫ��Ԫ
error2 = zeros(3,5);   % ��������

% 
for index = 1:5 
    n = num(index);
    
    A0 = zeros(n,n); b0 = zeros(n,1); % A0:ԭʼ�����еľ��� b0:ԭʼ�����е�����
    realSol = ones(n,1); % ��ȷ������
    % ����ԭʼ������������ֵ
    A0(1,1) = 6; A0(2,1) = 8; A0(n-1,n) = 1; A0(n,n) = 6; 
    b0(1,1) = 7; b0(n,1) = 14;
    for i = 2:n-1 
        A0(i-1,i) = 1;
        A0(i,i) = 6;
        A0(i+1,i) = 8;
        b0(i,1) = 15;
    end

    % 1 ֱ��Guass��ȥ��
    A = LU_norm(A0); b = b0; % LU_norm: ֱ�ӽ���Guass��ȥ��
    for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1) = b(1)/A(1,1);
        % ������
    error1(1,index) = norm(b - realSol);      % ������뾫ȷ������2���� 
    error2(1,index) = max(abs(b - realSol));  % ������뾫ȷ����������� 
    
    % 2 ����ԪGuass��ȥ��
    [A,u] = LU_col(A0); b = b0; % LU_col: ����ԪGuass��ȥ�� ����A=LU ����u(�����н�����Ϣ)
        
    for k = 1:n-1  % ���������н���( ����P^(-1)b )
        mid = b(k);
        b(k) = b(u(k));
        b(u(k)) = mid;
    end
    for i = 1:n-1  % ǰ������ Ly=b(y��¼��b��)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2 % �ش����� Ux=y(x��¼��b��)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1) = b(1)/A(1,1);
        % ������
    error1(2,index) = norm(b - realSol);     % ������뾫ȷ������2����
    error2(2,index) = max(abs(b - realSol)); % ������뾫ȷ�����������
    
    % 3 ȫ��ԪGuass��ȥ��
    [A,u,v] = LU_all(A0); b = b0; %LU_all: ȫ��ԪGuass��ȥ�� ����A=LU ����u(�����н�����Ϣ) ����v(�����н�����Ϣ)

    for k = 1:n-1    % ���������н���( ����P^(-1)b )
        mid = b(k);
        b(k) = b(u(k));
        b(u(k)) = mid;
    end
    for i = 1:n-1    % ǰ������ Ly=b(y��¼��b��)
        b(i+1:n) = b(i+1:n) - b(i)* A(i+1:n,i);
    end
    for j = n:-1:2   % �ش����� Ux=y(x��¼��b��)
        b(j) = b(j)/A(j,j);
        b(1:j-1) = b(1:j-1) - b(j)* A(1:j-1,j);
    end
    b(1,1) = b(1,1)./A(1,1);
    for k = n-1:-1:1 %���������н���(����Q^(-1)b )
        mid = b(k);
        b(k) = b(v(k));
        b(v(k)) = mid;
    end
        % �������
    error1(3,index) = norm(b - realSol);     % ������뾫ȷ������2����
    error2(3,index) = max(abs(b - realSol)); % ������뾫ȷ�����������
end
%% ������
error1
error2
