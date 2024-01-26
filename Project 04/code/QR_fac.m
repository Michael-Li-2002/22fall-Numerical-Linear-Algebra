function [A, d] = QR_fac(A)
% ����Householder�任��QR�ֽ�
% try: ��񻯴洢��������
[m,n] = size(A);
d = zeros(n,1); % �洢 Householder �任�ж�Ӧ��ϵ��
for j = 1:n
    if j < m    
        [v,beta] = house(A(j:m,j)); % ����Ӧ���е� Householder �任
        % ���Ӿ������ Householder �任
        % ������������㻯Ϊ��������������������������
        A(j:m,j:n) = A(j:m, j:n) - (beta * v) * (v' * A(j:m, j:n)); 
        d(j) = beta * v(1)^2;         % �洢��񻯺��ϵ�� beta
        A(j+1:m,j) = v(2:m-j+1)/v(1); % �洢��񻯺����������ķ���
    end
end

