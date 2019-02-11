%Make gamma likelihood heatmap equivalent maximizing mixture

sigma = 0.4;
phi = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));
Y = [1;2];

%Solution
Q = MixtureLikelihoodMovingMasses2(phi,Y);

options.save = true;
options.plot_Q = true;
options.plot_Q_unsimplified = false;
options.xx_padding = 1;
options.filename = 'gamma_likelihood_heatmap_solution.png';
options.main_colour = 'm';
options.pos = [200 200 600 450];

fig = plot_mixture(phi, Q, Y, options);
