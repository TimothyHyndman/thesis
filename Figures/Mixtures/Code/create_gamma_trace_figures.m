%Parameters and function defs
sigma = 0.4;
f = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));
Y = [1;2];
theta_res = 800;
george = linspace(-10,10,theta_res);
fred = f(Y,george);
x1max = max(fred(1,:));
x2max = max(fred(2,:));


main_colour = [0 51 104]/255;
marker_colour = [255 189 76]/255;
marker_colour = 'r';
marker_size = 80;

fig = figure('pos',[200, 200, 450, 450]);

theta_values = linspace(-10,1,theta_res);
gamma = f(Y,theta_values);
h = plot(gamma(1,:),gamma(2,:),'k','LineWidth',1.5');
xlabel('\gamma_1', 'FontSize', 20)
ylabel('\gamma_2', 'FontSize', 20)
axis equal
axis([0,x1max,0,x2max])
% export_fig GammaTrace0331.pdf


% set(gca, 'units', 'normalized'); %Just making sure it's normalized
% Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
%                                  %[Left Bottom Right Top] spacing
% NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
% set(gca, 'Position', NewPos);

saveas(h, 'GammaTrace01.png')

% pause

fig = figure('pos',[200, 200, 600, 450]);
xx = linspace(-1,4,100);
theta_loc = 1;
bob = f(xx',theta_loc);
h = plot(xx,bob, 'Color', main_colour, 'LineWidth',4);
hold on
scatter(Y,f(Y,theta_loc),marker_size, 'filled', 'MarkerFaceColor', marker_colour);
hold off

xlabel('x', 'FontSize', 20)
ylabel('f(x;\theta)', 'FontSize', 20)

% export_fig GammaTraceDensity01.pdf
saveas(h, 'GammaTraceDensity01.png')

% pause
fig = figure('pos',[200, 200, 450, 450]);
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

fig = figure('pos',[200, 200, 600, 450]);
xx = linspace(-1,4,100);
theta_loc = 2;
bob = f(xx',theta_loc);
h = plot(xx,bob, 'Color', main_colour, 'LineWidth',4);
hold on
scatter(Y,f(Y,theta_loc),marker_size,'filled', 'MarkerFaceColor', marker_colour);
hold off
% export_fig GammaTraceDensity02.pdf
xlabel('x', 'FontSize', 20)
ylabel('f(x;\theta)', 'FontSize', 20)
saveas(h, 'GammaTraceDensity02.png')

% pause
fig = figure('pos',[200, 200, 450, 450]);
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

fig = figure('pos',[200, 200, 600, 450]);
xx = linspace(-1,4,100);
theta_loc = 3;
bob = f(xx',theta_loc);
h = plot(xx,bob, 'Color', main_colour, 'LineWidth',4);
hold on
scatter(Y,f(Y,theta_loc),marker_size,'filled', 'MarkerFaceColor', marker_colour);
hold off
% export_fig GammaTraceDensity03.pdf
xlabel('x', 'FontSize', 20)
ylabel('f(x;\theta)', 'FontSize', 20)
saveas(h, 'GammaTraceDensity03.png')