clc;
close all;
clear;
format longEng

N = 1000; 
L = 1e-10; 
h = 6.62607015e-34; 
hbar = h / (2 * pi); 
m = 9.1093837e-31; 
ev = 1.6e-19; 
f0 = 2; 
K = m * (2 * pi * f0)^2 / ev; 

dx = L / N;
x = -2 * L * 1e10 : 1 / N : 2 * L * 1e10 - 1 / N; 
v = 0.5 * K * x.^2; 

H = zeros(4 * N, 4 * N);

op = -hbar^2 / (2 * m * dx^2) * [1, -2, 1];

for i = 2 : 4 * N - 1
    H(i, i) = op(2) + v(i);
    H(i, i - 1) = op(1);
    H(i, i + 1) = op(3);
end
H(1, 1) = op(2) + v(1);
H(1, 2) = op(3);
H(4 * N, 4 * N) = op(2) + v(4 * N);
H(4 * N, 4 * N - 1) = op(1);

[Psi, E] = eig(H);
E = diag(E);

figure
hold on
leg = [];
ran = [1, 25, 100];
for k = 1 : length(ran)
    [En, n] = mink(E, ran(k));
    plot(x, -Psi(:, n(end)) + 0.04 * ran(k))
    leg{k} = sprintf('\\Psi_{%d}', ran(k));
end
legend(leg)
title('Eigenfunctions')
xlabel('Position (x)')
ylabel('Wavefunction (\Psi)')

n = 1 : 100;
figure;
plot(mink(E, 100), 'o');
hold on
E_T = (n + 0.5) * hbar * 2 * pi * f0 / ev;
plot(E_T);
title('Energy Levels')
xlabel('n')
ylabel('Energy (E)')
legend('Numerical', 'Theoretical')

