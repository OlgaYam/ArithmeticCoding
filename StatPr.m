P = dlmread('matrixDistr.dat');
M = size(P,1);

% ���������� ������������ ������������
% ���������������� �������
B = P;
P = P';
% ��������� 1 �� ������������ ���������
for i = 1:M
    P(i,i) = P(i,i) - 1;
end
% �������� ��������� ������� ��������� (����� ���� ������������ = 1)
P(M,:) = ones(1,M);
% ������� �������� �������
P = inv(P);
% �������� �����������
Pr = P(:, M);

dlmwrite('Probabilities.txt',Pr);

H = 0;
h = zeros(1, M);
for i = 1:M
    for j = 1: M
        if(B(i,j) ~= 0) 
            h(j) = h(j) + log2(B(i,j))*B(i,j);
        end
    end
    H = H + h(i)*Pr(i);
end

H = H *(-1);
H