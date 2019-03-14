% m-point example

sigma = 1;
phi = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));
%since sigma = 1...
phidash = @(Y,theta) (1/(sqrt(2*pi)))*(repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])) .* exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./ 2);
%
phidashdash = @(Y,theta) (1/(sqrt(2*pi)))*((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2 - 1) .* exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./ 2);

rng(1)
X = [-2;2];
[Q2, likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X);

Q1.Support = 0;
Q1.ProbWeights = 1;

f_Q = Q1.ProbWeights*phi(X,Q1.Support)';

theta = 0;
fdashdash = phidashdash(X, theta)';

Ddashdash = sum(fdashdash ./ f_Q, 2)

fdash = phidash(X, theta)';

fred = sum((fdash ./ f_Q).^2,2)

f = phi(X, theta)'
george = sum((f ./ f_Q).^2, 2)



DQ = @(theta) D(phi, X, theta, Q1);


thetas = linspace(0, 2, 1000);
ds = DQ(thetas);


%numerically confirming ddashdash

dtheta = thetas(2) - thetas(1);
ds1 = (ds(2) - ds(1)) / dtheta

ds2 = (ds(3) - ds(2)) / dtheta

Ddd = (ds2 - ds1) / dtheta


plot(thetas, ds)

function derivative = D(phi,X,theta,Q)
f_Q = Q.ProbWeights*phi(X,Q.Support)';
f_theta = phi(X,theta)';
denom = repmat(f_Q,[length(theta),1]);
derivative = sum(f_theta./denom - 1,2);
end

