% Fixed masses example

%Initiate random seeds to be able to replicate the results
rng(5000)

% X
n = 500;
dchi = 3;
X = chi2rnd(dchi,1,n);
varX = 2 * dchi;
X = X / sqrt(varX);
varX = 1;
varXlong=2*dchi;
cvar=sqrt(varXlong);
truedens=@(xx) cvar*chi2pdf(xx*cvar,dchi);

% U
NSR = 1;
sigU = sqrt(NSR * varX);
U = normrnd(0, sigU, 1, n);

% W
W = X + U;
h=bwsjpiSM(W');


xx_fixed = linspace(min(W), max(W), 100);

[yy_fixed, Q_fixed, tt_fixed, optim_values_fixed, misc_variables_fixed] = fixed_mass_deconvolve(W, xx_fixed);

options.save = false;
options.filename = 'fixed_masses_example_NSR1.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = true;
options.naivebw = h;
options.plot_naive = true;
plot_deconvolution_graph(Q_fixed, xx_fixed, yy_fixed, W, options);