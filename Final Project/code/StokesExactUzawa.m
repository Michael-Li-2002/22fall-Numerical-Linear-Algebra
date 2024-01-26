%% ��ֵ��������ҵ �ڶ���(Exact Uzawa�������Stokes���̣�
% ���ò������������ʵĲ����󵥻����м���
%% ���ò���
level = 7; alpha = 1; % �����ģ N=2^level; ���������� alpha
%% ���÷����鳣���������
N = 2^level; h = 1/N; % ����λ���� h
F = zeros(N,N-1); G = zeros(N-1,N); D = zeros(N);
for j=1:N-1 % �������Է�����ĳ�����
    F(1,j) = f(j*h,0.5*h) - 2*pi*(1-cos(2*pi*j*h))*N;
    F(N,j) = f(j*h,1-0.5*h) + 2*pi*(1-cos(2*pi*j*h))*N;
    G(j,1) = g(0.5*h,j*h) + 2*pi*(1-cos(2*pi*j*h))*N;
    G(j,N) = g(1-0.5*h,j*h) - 2*pi*(1-cos(2*pi*j*h))*N;
    for i=2:N-1
        F(i,j) = f(j*h,(i-0.5)*h);
        G(j,i) = g((i-0.5)*h,j*h);
    end
end
real_U = zeros(N,N-1); real_V = zeros(N-1,N); real_P = zeros(N); 
for i=1:N % �������
    for j = 1: N-1
        real_U(i,j) = (1-cos(2*pi*j*h))*sin(2*pi*(i-0.5)*h);
        real_V(j,i) = -(1-cos(2*pi*j*h))*sin(2*pi*(i-0.5)*h);
    end
    for j = 1:N
        real_P(i,j) = (j-0.5)^3*h^3/3-1/12;
    end
end
%% ���Stokes����
U = zeros(N,N-1); V = zeros(N-1,N); P = zeros(N);
[F_err,G_err,D_err] = residue(U,V,P,F,G,D); % ���ʼ����
r0 = sqrt(sum(sum(F_err.^2))+sum(sum(G_err.^2))+sum(sum(D_err.^2))); % ��ʼ����2���� r0
res = r0; cycle = 0; % ��ǰ����2���� res; ��ѭ������ cycle
while res/r0 > 10^(-8)
    cycle = cycle + 1;
    [U,V] = CGforStokes(U,V,P,F,G); % �����ݶȷ����������
    [~,~,D_err] = residue(U,V,P,F,G,D); % ����B^TU
    P = P + alpha * D_err; % ����ѹ��
    [F_err,G_err,D_err] = residue(U,V,P,F,G,D); % ���㵱ǰ����
    res = sqrt(sum(sum(F_err.^2))+sum(sum(G_err.^2))+sum(sum(D_err.^2)));
    disp(['ExactUzawa complete ',num2str(cycle),' times, with residue ',num2str(res)])
end
%% ������
error = h*sqrt(sum(sum((real_U-U).^2)) + sum(sum((real_V-V).^2)));
disp(['Numerical error: ',num2str(error),'.'])