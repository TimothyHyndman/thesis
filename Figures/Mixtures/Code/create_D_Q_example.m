
function create_D_Q_example()
%Parameters and function defs

sigma = 0.5;
phi = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));

rng(1)
X = unifrnd(-3,3,[50,1]);
[Q, likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X);




DQ = @(theta) D(phi,X,theta,Q);
res = 1000;
xx = linspace(min(X),max(X),res);
yy = DQ(xx);

figure('pos', [200 200 800 450])
fig = plot(xx,yy, 'Color','k','LineWidth',3);
xlabel('\theta','FontSize',20)

saveas(fig, 'D_Q_example.png')


options.plot_Q = true;
options.plot_Q_unsimplified = false;
options.save = true;
options.filename = 'D_Q_example_mixture.png';
options.xx_padding = 1;

plot_mixture(phi, Q, X, options, Q_unsimplified);


end

function derivative = D(phi,X,theta,Q)
f_Q = Q.ProbWeights*phi(X,Q.Support)';
f_theta = phi(X,theta)';
denom = repmat(f_Q,[length(theta),1]);
derivative = sum(f_theta./denom - 1,2);
end