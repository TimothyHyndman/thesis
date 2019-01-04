% Fixed masses example

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


xx_emp = linspace(min(W), max(W), 100);

Q_emp.Support = W;
Q_emp.ProbWeights = ones(1, n)/n;

length_tt = 100;
mu_K2 = 6;
RK = 1024 / 3003 / pi;
hnaive = ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(var(W)) * n^(-1/5);
hmin = hnaive/3;
tt_emp = linspace(-1/hmin, 1/hmin, length_tt);

[rehatphiW, imhatphiW] = compute_phi_W(tt_emp, W);
normhatphiW = sqrt(rehatphiW.^2 + imhatphiW.^2)';
t_star = find_t_cutoff(normhatphiW, tt_emp, n);

tt_emp = linspace(-t_star, t_star, length_tt);

[rehatphiW, imhatphiW] = compute_phi_W(tt_emp, W);
hat_phi_W = complex(rehatphiW, imhatphiW).';
[~, ~, sqrt_psi_hat_W] = compute_psi_W(tt_emp, W);
sqrt_psi_hat_W = sqrt_psi_hat_W';

weight = KernelWeight('Epanechnikov',tt_emp);

tp_emp = calculate_tp(tt_emp, Q_emp.ProbWeights, Q_emp.Support, hat_phi_W, sqrt_psi_hat_W, weight);
    

yy_emp = decon_err_sym_pmf2pdf(xx_emp, tt_emp, Q_emp.Support, Q_emp.ProbWeights, W, []);

options.save = false;
options.filename = 'emp_masses_example.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = true;
options.naivebw = h;
options.plot_naive = true;
plot_deconvolution_graph(Q_emp, xx_emp, yy_emp, W, options);