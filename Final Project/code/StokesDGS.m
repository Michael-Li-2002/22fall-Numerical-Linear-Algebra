%% ��ֵ��������ҵ ��һ��(��DGSΪĥ���ӵĶ������񷽷����Stokes���̣�
% ���ò������������ʵĲ����󵥻����м���
%% ���ò���
level = 11; r = 10;   % �����ģ N=2^level; �������Ϊ r+1(L=2^r)
nu1 = 2; nu2 = 2;   % ǰĥ����� nu1; ��ĥ����� nu2
smoother = @DGS;    % ĥ���� DGS
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
for i=1:N  % �������
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
r0 = 0;
for i=1:N
    for j=1:N-1
       r0 = r0 + F(i,j)^2 + G(j,i)^2;
    end
end
r0 = sqrt(r0); res = r0; cycle = 0; % ��ʼ����2���� r0; ��ǰ����2���� res; ѭ������ cycle
while res/r0 > 10^(-8) 
    cycle = cycle + 1;
    [U,V,P] = Vcycle(level,r,nu1,nu2,U,V,P,F,G,D,smoother); % ����һ��V-cycle
    [F_err,G_err,D_err] = residue(U,V,P,F,G,D); % �������
    res = 0;
    for i = 1:N
        for j = 1:N-1
            res = res + F_err(i,j)^2 + G_err(j,i)^2;
        end
        for j = 1:N
            res = res + D_err(i,j)^2;
        end
    end
    res = sqrt(res); % �������2����
    disp(['vcycle complete ',num2str(cycle),' times, with residue ',num2str(res)])
end
%% ������
error = h*sqrt(sum(sum((real_U-U).^2)) + sum(sum((real_V-V).^2)));
disp(['Numerical error: ',num2str(error),'.'])