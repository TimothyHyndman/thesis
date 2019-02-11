%Normal
rng(8)
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));

X = unifrnd(-2,2,[50,1]);

[Q, likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X);

options.plot_Q = true;
options.plot_Q_unsimplified = true;
options.save = false;
options.filename = 'collapsing_distribution.png';

plot_mixture(phi, Q, X, options, Q_unsimplified);