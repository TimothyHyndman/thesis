rng(4)
%Normal
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));
%Uniform
% sigma = 1;
% phi = @(x,u) (1/2*sigma)*(abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))<=sigma);
%Triangle
% sigma = 1;
% phi = @(x,u) (1/sigma)*(1 - (1/sigma)*abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))).*(abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))<=sigma);
%Cauchy
% lambda = sqrt(3);   %makes inflection points at +- 1
% phi = @(x,u) (1/(pi*lambda))*(1 + ((repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))./lambda).^2).^(-1);

X = chi2rnd(3, [500,1]);

[Q, likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X);

options.plot_Q = true;
options.plot_Q_unsimplified = false;
options.save = true;
options.filename = 'chi2_n500_motivation.png';

plot_mixture(phi, Q, X, options, Q_unsimplified);