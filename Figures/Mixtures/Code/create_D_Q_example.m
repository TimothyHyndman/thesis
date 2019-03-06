%Parameters and function defs
sigma = 1;
phi = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));

rng(1)
X = unifrnd(-4,4,[500,1]);
[Q, likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X);

options.plot_Q = true;
options.plot_Q_unsimplified = true;
options.save = false;
options.filename = 'blahblahblah.png';

plot_mixture(phi, Q, X, options, Q_unsimplified);