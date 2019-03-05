sigma = 1;
% Normal Density
% g = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));

% f = @(Y, theta) 0.4*g(Y, theta - 4) + 0.4*g(Y, theta + 4) + 0.2*g(Y, theta);
% Y = [0,3]';

% Custom
f = @(Y, theta) max(0, 1 - 2*abs(repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])));
Y = [0;0.4];

% f = @(Y, theta) (abs(repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])) <= 1) .* (1 + cos(3 * pi * ( repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1]) ) ));

xx = linspace(-1,1, 1000)';
yy = f(xx, 0);
figure(1)
plot(xx, yy);

%Cauchy density
% lambda = sqrt(3);   %makes inflection points at +- 1
% f = @(x,u) (1/(pi*lambda))*(1 + ((repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))./lambda).^2).^(-1);

%Gamma
theta_res = 8000;
theta_values = linspace(-10,10,theta_res);
% Y = [0;0.75];
gamma = f(Y,theta_values);
% k = convhull(gamma(1,:),gamma(2,:));

% %Optimization surface
% x1max = max(gamma(1,:));
% x2max = max(gamma(2,:));
% xres = 1000;
% yres = 1000;
% x = linspace(0,x1max,xres);
% y = linspace(0,x2max,yres);
% for i = 1:xres
%     for j = 1:yres
%         S(i,j) = x(i)*y(j);
%     end
% end
% imagesc([0,x1max],[0,x2max], S);
% set(gca,'YDir','normal')

%Optimal solution and figure 1
plotflag = 1;
% [solution_grid_points,solution_masses] = MixtureLikelihoodMovingMasses(f,plotflag,Y);
Q = MixtureLikelihoodMovingMasses2(f,Y);

% options.xx_padding = 1;
% options.plot_data_points = true;
% options.filename = 'triangle_density_mixture.png';
% options.save = false;
% options.plot_Q = true;
% plot_mixture(f, Q, Y, options)

optimal_thetas = f(Y,Q.Support);
optimal_u = Q.ProbWeights*optimal_thetas';

a = optimal_u(1)*2;
p = 0.92;
hyperplane_x = linspace((1-p)*a, p*a, 2);
hyperplane_y = linspace(p*a, (1-p)*a, 2);

% export_fig MixingSolSigma06.pdf

%Plot
h = figure(2);
axis equal tight
clf
hold on
% contourf(x,y,S,200,'LineStyle','none')
% imagesc([0,x1max],[0,x2max], S);
set(gca,'YDir','normal')
plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
plot(hyperplane_x, hyperplane_y,'--','LineWidth',1.5');
% plot(gamma(1,k),gamma(2,k),'k','LineStyle','--','LineWidth',1.5');
% scatter(optimal_thetas(1,:),optimal_thetas(2,:),'m','filled')
scatter(optimal_u(1),optimal_u(2),'b','filled')
hold off
xlabel('\gamma_1','FontSize',20)
ylabel('\gamma_2','FontSize',20)
axis equal
axis([0, a+0.01, 0, a])
% export_fig GammaSigma06.pdf

filename = 'triangle_density_gamma.png';
set(h, 'PaperPositionMode', 'auto');
saveas(h, filename)