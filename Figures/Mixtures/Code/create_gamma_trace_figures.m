%Parameters and function defs
sigma = 0.4;
f = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));
Y = [1;2];
theta_res = 800;
george = linspace(-10,10,theta_res);
fred = f(Y,george);
x1max = max(fred(1,:));
x2max = max(fred(2,:));


%Start plotting
% theta_values = linspace(-10,0,theta_res);
% gamma = f(Y,theta_values);
% plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
% xlabel('u_1')
%     ylabel('u_2')
%     axis equal
%     axis([0,x1max,0,x2max])
% % export_fig GammaTrace00.pdf
% 
% pause
% 
% xx = linspace(-1,4,100);
% theta_loc = 0;
% bob = f(xx',theta_loc);
% plot(xx,bob)
% hold on
% scatter(Y,f(Y,theta_loc));
% hold off
% % export_fig GammaTraceDensity00.pdf

% pause

theta_values = linspace(-10,1,theta_res);
gamma = f(Y,theta_values);
h = plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
xlabel('\gamma_1', 'FontSize', 20)
ylabel('\gamma_2', 'FontSize', 20)
axis equal
axis([0,x1max,0,x2max])
export_fig GammaTrace0331.pdf
% saveas(h, 'GammaTrace01.png')

% pause

xx = linspace(-1,4,100);
theta_loc = 1;
bob = f(xx',theta_loc);
h = plot(xx,bob);
hold on
scatter(Y,f(Y,theta_loc));
hold off
% export_fig GammaTraceDensity01.pdf
saveas(h, 'GammaTraceDensity01.png')

% pause

theta_values = linspace(-10,2,theta_res);
gamma = f(Y,theta_values);
h = plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
xlabel('\gamma_1', 'FontSize', 20)
ylabel('\gamma_2', 'FontSize', 20)
axis equal
axis([0,x1max,0,x2max])
% export_fig GammaTrace02.pdf
saveas(h, 'GammaTrace02.png')

% pause


xx = linspace(-1,4,100);
theta_loc = 2;
bob = f(xx',theta_loc);
h = plot(xx,bob);
hold on
scatter(Y,f(Y,theta_loc));
hold off
% export_fig GammaTraceDensity02.pdf
saveas(h, 'GammaTraceDensity02.png')

% pause

theta_values = linspace(-10,10,theta_res);
gamma = f(Y,theta_values);
h = plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
xlabel('\gamma_1', 'FontSize', 20)
ylabel('\gamma_2', 'FontSize', 20)
axis equal
axis([0,x1max,0,x2max])
% export_fig GammaTrace03.pdf
saveas(h, 'GammaTrace03.png')
% pause

xx = linspace(-1,4,100);
theta_loc = 3;
bob = f(xx',theta_loc);
h = plot(xx,bob);
hold on
scatter(Y,f(Y,theta_loc));
hold off
% export_fig GammaTraceDensity03.pdf
saveas(h, 'GammaTraceDensity03.png')