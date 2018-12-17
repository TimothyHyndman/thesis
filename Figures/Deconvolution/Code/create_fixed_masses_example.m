% Fixed masses example

%Initiate random seeds to be able to replicate the results
% myseed=5000;
% rand('state',myseed);
% randn('state',myseed);

% X
NSR = 0.2;
n = 500;
dchi = 3;
X = chi2rnd(dchi,1,n);
varX = 2 * dchi;
X = X / sqrt(varX);
varX = 1;

% U
sigU = sqrt(NSR * varX);
U = normrnd(0, sigU, 1, n);

W = X + U;
xx = linspace(min(W), max(W), 100);

[yy, Q] = fixed_mass_deconvolve(W, xx);