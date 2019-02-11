
sigma = 0.4;
% Normal Density
f = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));

%Gamma
theta_res = 800;
theta_values = linspace(-10,10,theta_res);
Y = [1;2];
gamma = f(Y,theta_values);
k = convhull(gamma(1,:),gamma(2,:));

%Optimization surface
x1max = max(gamma(1,:));
x2max = max(gamma(2,:));
xres = 1000;
yres = 1000;
x = linspace(0,x1max,xres);
y = linspace(0,x2max,yres);
for i = 1:xres
    for j = 1:yres
        S(i,j) = x(i)*y(j);
    end
end
imagesc([0,x1max],[0,x2max], S);
set(gca,'YDir','normal')

%Optimal solution and figure 1
plotflag = 1;
% [solution_grid_points,solution_masses] = MixtureLikelihoodMovingMasses(f,plotflag,Y);
Q = MixtureLikelihoodMovingMasses2(f,Y);
optimal_thetas = f(Y,Q.Support);
optimal_u = Q.ProbWeights*optimal_thetas';

% export_fig MixingSolSigma06.pdf

%Plot
h = figure(2);
axis equal tight
clf
hold on
% contourf(x,y,S,200,'LineStyle','none')
imagesc([0,x1max],[0,x2max], S);
set(gca,'YDir','normal')
plot(gamma(1,:),gamma(2,:),'w','LineWidth',1.5');
plot(gamma(1,k),gamma(2,k),'k','LineStyle','--','LineWidth',1.5');
scatter(optimal_thetas(1,:),optimal_thetas(2,:),'m','filled')
scatter(optimal_u(1),optimal_u(2),'y','filled')
hold off
xlabel('\gamma_1','FontSize',20)
ylabel('\gamma_2','FontSize',20)
axis equal
axis([0,x1max,0,x2max])
% export_fig GammaSigma06.pdf

filename = 'gamma_likelihood_heatmap.png';
% saveas(h, filename)