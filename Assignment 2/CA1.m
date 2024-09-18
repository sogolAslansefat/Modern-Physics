clc
clear
close all; 
format long
eV = 1.602e-19;
L = 1e-10;
me = 9.11e-31;
hbar=1.055e-34;
const = -(hbar)^2/(2*me);
U = @(x) 0*x;
deltaX = 0.001*1e-10;
x = (0 : deltaX : 1e-10)';
[n, ~] = size(x);
ham = (const/(deltaX^2))*ones(n-1, 1)*[1, -2, 1];
ham(:,1) = ham(:,1) + U(x(1:end - 1));
part1 = toeplitz([1, zeros(1,n-2)], [1, zeros(1,n)]);
part2 = toeplitz(zeros(1,n-1), [0, 1, zeros(1,n-1)]) ;
part3 = toeplitz(zeros(1,n-1), [0, 0, 1, zeros(1,n-2)]) ;
B = (part1.*ham(:, 1));
C = (part2.*ham(:, 2));
D = (part3.*ham(:, 3));
A = B + C + D;
HamMat = A(:,2:end - 1);
[psi,En] = eig(HamMat);
[~, nn] = size(En);
En = (En*ones(nn, 1))/eV;
K = @(U, E) sqrt(2*me*(E-U)/(hbar)^2);
kn = K(U(x(2:end)), En*eV);
nL=1:100;
scatter(nL, En(1:100)*eV)
hold on
E_th = nL.^2*pi^2*hbar^2 / (2*me*L^2);
plot(nL,E_th);
figure
plot(x(1:end-1), -psi(:, 1))
hold on
plot(x(1:end-1), -psi(:, 25))
leg{1} = sprintf('\\Psi_{%d}', 1);
leg{2} = sprintf('\\Psi_{%d}', 25);
legend(leg)








