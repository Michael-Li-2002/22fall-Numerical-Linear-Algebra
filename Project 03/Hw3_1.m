%% ��ֵ���� ��3���ϻ���ҵ ��� 2100010793
% P1: ���� 5-20 �� Hilbert �����������
% ֱ�ӵ�������С� �ڹ������鿴����'cond_num'�м��������������ֵ
%% ���þ���
cond_num = zeros(1,16);
%% ����������
for n = 5:20
    A = hilb(n);
    InfNormA = sum(A(1,:));
    InfNormAInv = InftyNorm_Inv(A);
    cond_num(n-4) = InfNormA * InfNormAInv;
end