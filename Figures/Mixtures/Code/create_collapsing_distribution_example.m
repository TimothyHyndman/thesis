%Normal
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));

X = unifrnd(-2,2,[50,1]);

[Q,likelihood] = MixtureLikelihoodMovingMasses2(phi,X);

scatter(Q.Support, Q.ProbWeights)