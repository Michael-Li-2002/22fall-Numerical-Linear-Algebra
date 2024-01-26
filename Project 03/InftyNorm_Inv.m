function [norm1] = InftyNorm_Inv(A)
% �������A �������A^-1�������
OK = 1; % �����ж��Ƿ���ֹѭ��
[~,n] = size(A);
[A0,u] = LU_col(A);% ��ǰ������Ԫ���Ƿֽ�
x = ones(n,1)/n; B0 = A0'; % �����ʼ���� �Ծ�����ת��
while OK == 1 
    w = x;  
    % ��� A^Tw=x
    for i = 1:n-1  % ǰ������ U^Ty=b(y��¼��b��)
        w(i) = w(i)/B0(i,i);
        w(i+1:n) = w(i+1:n) - w(i)* B0(i+1:n,i);
    end
    w(n) = w(n)/B0(n,n);
    for j = n:-1:2 % �ش����� L^Tx=y(x��¼��b��)
        w(1:j-1) = w(1:j-1) - w(j)* B0(1:j-1,j);
    end
    for k = n-1:-1:1  % ���������н���( ����P^Tb )
        mid = w(k);
        w(k) = w(u(k));
        w(u(k)) = mid;
    end
    z = sign(w);
    % ��� Az=v
    for k = 1:n-1  % ���������н���( ����P^(-1)b )
        mid = z(k);
        z(k) = z(u(k));
        z(u(k)) = mid;
    end
    for i = 1:n-1  % ǰ������ Ly=z(y��¼��z��)
        z(i+1:n) = z(i+1:n) - z(i)* A0(i+1:n,i);
    end
    for j = n:-1:2 % �ش����� Ux=y(x��¼��z��)
        z(j) = z(j)/A0(j,j);
        z(1:j-1) = z(1:j-1) - z(j)* A0(1:j-1,j);
    end
    z(1) = z(1)/A0(1,1);
    inftynorm = 0; % z�������; 
    pos = 1; % zĳ��������ģ���������ȵ�����
    for index = 1:n
        if abs(z(index)) > inftynorm
            inftynorm = abs(z(index));
            pos = index;
        end
    end
    if inftynorm > z'* x % �������ֵ ��������
        x = zeros(n,1);
        x(pos) = 1;
    else % ֹͣ����
        norm1 = sum(abs(w));
        OK = 0;
    end
end
    

