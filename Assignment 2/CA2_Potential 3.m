clc;
close all;
clear;
format longEng;

N = 1000;
L = 1e-10;
h = 6.62607015e-34;
hbar = h / (2 * pi);
m = 9.11e-31;
ev = 1.6e-19;
V = 1e9 * ev;

dx = L / N;
x = linspace(-L, 2*L, 3*N);
v_ramp = (x(N+1:end) + L) * V / (2 * L);

H = zeros(3 * N);

op = -hbar^2 / (2 * m * dx^2) * [1, -2, 1];
Vmat = V * eye(3 * N);

% ????? ??????? ??? ??? ?? ????? ??? ? ???
Vmat(N + 1:3 * N, N + 1:3 * N) = diag(v_ramp);

% ?? ???? ?????? ???????? ?? ??????? ????????
for i = 1:3 * N
    if i > 1
        H(i, i - 1) = op(1);
    end
    H(i, i) = op(2);
    if i < 3 * N
        H(i, i + 1) = op(3);
    end
end

Hamiltonian = H + Vmat;

[Psi, E] = eig(Hamiltonian);
E = diag(E);

figure;
hold on;
leg = {};
ran = [1, 25];
for k = 1:length(ran)
    [En, n] = mink(E, ran(k));
    plot(x, Psi(:, n(end)) + 0.04 * ran(k));
    leg{k} = sprintf('\\Psi_{%d}', ran(k));
end
legend(leg);
xlabel('x');
ylabel('Wavefunction');
title('Wavefunctions');

n = 1:100;
figure;
hold on;
plot(mink(E, 100), '--');
xlabel('Index');
ylabel('Energy (J)');
title('Energy levels');