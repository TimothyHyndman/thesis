% Moving masses example

%Initiate random seeds to be able to replicate the results
seed = 5000;
NSR = 0.2;
n = 5000;
dist_type = 'discrete2';
error_type = 'norm';

[W, truedens, X, U, truepmf] = generatedata(n,NSR,dist_type,error_type,seed);
h=bwsjpiSM(W');


xx_moving = linspace(min(W), max(W), 100);

decon_options.penalties = true;
[yy_moving, Q_moving, tt_moving, optim_values_moving] = moving_mass_deconvolve(W, xx_moving, 20, decon_options);

options.save = false;
options.filename = 'moving_masses_discrete_n5000_example.png';
options.plot_histogram = false;
options.plot_density = false;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = false;
options.true_pmf = truepmf;
options.plot_true_pmf = true;
options.naivebw = h;
options.plot_naive = false;
plot_deconvolution_graph(Q_moving, xx_moving, yy_moving, W, options);