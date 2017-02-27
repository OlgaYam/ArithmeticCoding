P = dlmread('matrixDistr.dat');
M = size(P,1);

% нахождение стационарных вероятностей
% транспонирование матрицы
B = P;
P = P';
% вычитание 1 из диагональных элементов
for i = 1:M
    P(i,i) = P(i,i) - 1;
end
% заменяем последнюю строчку единицами (сумма всех вероятностей = 1)
P(M,:) = ones(1,M);
% находим обратную матрицу
P = inv(P);
% получаем вероятности
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