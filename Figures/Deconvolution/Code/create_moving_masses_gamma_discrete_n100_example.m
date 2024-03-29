% Moving masses example

%Initiate random seeds to be able to replicate the results
seed = 5000;
NSR = 0.2;
n = 100;
dist_type = 'gamma';
error_type = 'discrete';

[W, truedens, X, U, truepmf] = generatedata(n,NSR,dist_type,error_type,seed);
h=bwsjpiSM(W');


xx_moving = linspace(min(W), max(W), 100);

decon_options.penalties = true;
[yy_moving, Q_moving, tt_moving, optim_values_moving] = moving_mass_deconvolve(W, xx_moving, 20, decon_options);

options.save = false;
options.filename = 'moving_masses_gamma_discrete_n100_example.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = true;
options.true_pmf = truepmf;
options.plot_true_pmf = false;
options.naivebw = h;
options.plot_naive = true;
plot_deconvolution_graph(Q_moving, xx_moving, yy_moving, W, options);