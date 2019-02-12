x2 = [1,2,3];

sigma = 1;
% Normal Density
f = @(Y,theta) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(Y,[1,length(theta)]) - repmat(theta,[length(Y),1])).^2)./(2*sigma^2));

%Gamma
theta_res = 800;
theta_values = linspace(-10,10,theta_res);

for l = 1:3
    Y = [0; x2(l)];
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
    Q = MixtureLikelihoodMovingMasses2(f,Y);
    optimal_thetas = f(Y,Q.Support);
    optimal_u = Q.ProbWeights*optimal_thetas';

    %Plot
    h = figure('pos',[200, 200, 450, 450]);
%     h = figure(2);
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

    filename = ['Sigma1x1_0-x2_',num2str(x2(l)),'.png'];
    
    set(gca, 'units', 'normalized'); %Just making sure it's normalized
    Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                     %[Left Bottom Right Top] spacing
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
    saveas(h, filename)
end