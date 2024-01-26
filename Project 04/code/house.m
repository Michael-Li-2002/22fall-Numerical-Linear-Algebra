function [v,beta] = house(x)
% ���� Householder �任 ��� H=I-beta*v*v^T �е� beta, v
% try: ��񻯴洢���� QR_fac ������
n = length(x);
x = x/max(abs(x));      % ��һ�������ֹ��������
sigma = x(2:n)'*x(2:n); 
v = x;
if sigma == 0    % ����������һ�������Ϊ 0 ��������� Householder �任
    beta = 0;    % �� beta Ϊ 0 ������ Householder �任
else
    alpha = sqrt(x(1)^2 + sigma); % ���� x �� ģ��
    if x(1) <= 0 % ����һ���������� 0 �����һ��ļ���
        v(1) = x(1) - alpha;
    else         % ����һ�������� 0 ��������㷽ʽ
        v(1) = -sigma/(x(1)+alpha);
    end          
    beta = 2 /(sigma + v(1)^2); % ����ϵ�� beta
end