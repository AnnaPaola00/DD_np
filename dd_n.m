clear all;
close all;
pkg load bim

L = 1.5e-6;
N = 100;
x = linspace(0,L,N)';
VA = -0.8;

D = 1e22 .* (x < .5*L);
A = 1e22 .* (x >= .5*L);
ni = 1e16;
q = 1.6e-19;
Vth = 26e-3;
epsilon = 4*8.85e-12;

% Guess iniziale
n1 = D(1)/2 + sqrt(D(1)^2 + 4*ni^2)/2;
p1 = ni^2/n1;
p2 = A(end)/2 + sqrt(A(end)^2 + 4*ni^2)/2;
n2 = ni^2/p2;

n = n1 .* (x <= .5*L) + n2 .* (x >  .5*L);
p = p2 .* (x >  .5*L) + p1 .* (x <= .5*L);

n = n1 .* (x <= .45*L) + n2 .* (x >  .45*L);
p = p2 .* (x >  .55*L) + p1 .* (x <= .55*L);


mup = 1e-1;
mun = .3e-1;

phi = Vth * log (n / ni);
phi(x>=.55*L) += VA;
phi0 = phi;

h = diff(x);

% Discretizzazione operatore di Laplace
P = bim1a_laplacian(x,epsilon,1);

% Matrice di Massa
M = bim1a_reaction(x,1,1);



[bp, bm] = bimu_bernoulli (diff (phi) / Vth);

% Discretizzazione equazione di continuità per gli elettroni

Cl = [-mun*Vth*(1./h).*bm; 0];
Cr = [0; -mun*Vth*(1./h).*bp];
Cd = mun*Vth*([(1./h).*bm; 0] + [0; (1./h).*bp]);
Cn = spdiags([Cl, Cd, Cr], -1:1, N,N);


Rp = P*phi + M*n;
Rc = Cn(2:end,2:end)*n(2:end);


