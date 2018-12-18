% Moving masses example

%Initiate random seeds to be able to replicate the results
rng(5000)

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
xx_moving = linspace(min(W), max(W), 100);

[yy_moving, Q_moving, tt_moving, optim_values_moving] = moving_mass_deconvolve(W, xx_moving);

options.save = false;
options.filename = 'moving_masses_example.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
plot_deconvolution_graph(Q_moving, xx_moving, yy_moving, W, options);