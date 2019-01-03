% VEM example

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
varXlong=2*dchi;
cvar=sqrt(varXlong);
truedens=@(xx) cvar*chi2pdf(xx*cvar,dchi);

% U
sigU = sqrt(NSR * varX);
U = normrnd(0, sigU, 1, n);

% W
W = X + U;
h=bwsjpiSM(W');


xx_VEM = linspace(min(W), max(W), 100);

[Q_VEM, tt_VEM, optim_values_VEM] = VEM_deconvolve(W, xx);

yy_VEM = decon_err_sym_pmf2pdf(xx_VEM, tt_VEM, Q_VEM.Support, Q_VEM.ProbWeights, W, []);

options.save = false;
options.filename = 'VEM_example.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = true;
options.naivebw = h;
options.plot_naive = true;
plot_deconvolution_graph(Q_VEM, xx_VEM, yy_VEM, W, options);
