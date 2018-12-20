% VEM like algorithm for TP

%Initiate random seeds to be able to replicate the results
rng(5000)

% X
NSR = 0.2;
n = 500;
dchi = 3;
X = chi2rnd(dchi,1,n);
varX = 2 * dchi;
X = X / sqrt(varX);
varX = 1;
varXlong=2*dchi;
cvar=sqrt(varXlong);
truedens=@(xx) cvar*chi2pdf(xx*cvar,dchi);

% U
sigU = sqrt(NSR * varX);
U = normrnd(0, sigU, 1, n);

% W
W = X + U;
h=bwsjpiSM(W');


xx = linspace(min(W), max(W), 100);

% ------------------------------------------------------------------------
varW = var(W);
dx=xx(2)-xx(1);
longt = 99;
mu_K2 = 6;
RK = 1024 / 3003 / pi;
hnaive = ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(var(W)) * n^(-1/5);
hmin = hnaive/3;
a = -1/hmin;
b = 1/hmin;
tt = a:((b - a) / longt):b;

% Compute phiW and psiW
[tt1, tt2, rehatphiW, imhatphiW, normhatphiW] = computephiW(tt,longt,W,n);
tt = tt1:(tt2-tt1)/longt:tt2;

[~, ~, sqrt_psi_hat_W] = computepsiW(tt, W, n);
hat_phi_W = complex(rehatphiW, imhatphiW);
weight = KernelWeight('Epanechnikov', tt);

% -------------------------------------------------------------------------
% Initial distribution
Q.Support = mean(W);
Q.ProbWeights = 1;

calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight)

looping = true;
counter = 0;
while looping
    % Add theta which minimizes tp to support
    theta_min = find_theta_star(Q, W, tt, hat_phi_W, sqrt_psi_hat_W, weight);
    Q.Support = [Q.Support, theta_min];
    Q.ProbWeights = [Q.ProbWeights, 0];
    
    %Rank support points by their tp_derivative value
    theta_derivatives = tp_derivative(Q.Support, Q, tt, hat_phi_W, sqrt_psi_hat_W, weight);
    [~, I] = sort(theta_derivatives);
    
    %Move mass from highest tp_masses to lowest
    objective_func = @(pj) calculate_tp(tt, pj, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight);
    Q = exchange_mass(Q, I(1), I, objective_func);
    
    % ALTERNATE BETWEEN RANKING BY DERIVATIVE AND RANKING BY ABSOLUTE SIZE
    % (to get rid of clusters and thin down the distribution a bit)
    
%     % Add theta which minimizes tp to support
%     theta_min = find_theta_star(Q, W, tt, hat_phi_W, sqrt_psi_hat_W, weight);
%     Q.Support = [Q.Support, theta_min];
%     Q.ProbWeights = [Q.ProbWeights, 0];
%     
%     %Rank point masses by size
%     [~, I] = sort(Q.ProbWeights, 'descend'); %Move smallest into largest
%     
%     %Move mass from smallest masses to theta_min
%     objective_func = @(pj) calculate_tp(tt, pj, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight);
%     Q = exchange_mass(Q, length(Q.ProbWeights), I, objective_func);
    
    counter = counter + 1;
    if counter > 100
        looping = false;
    end
end

calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight)

xx = linspace(min(W), max(W), 100);
yy = decon_err_sym_pmf2pdf(xx, tt, Q.Support, Q.ProbWeights, W, []);

options.save = false;
options.filename = 'VEM_example.png';
options.plot_histogram = false;
options.plot_density = true;
options.plot_masses = true;
options.true_dens = truedens;
options.plot_true_dens = true;
options.naivebw = h;
options.plot_naive = true;
plot_deconvolution_graph(Q, xx, yy, W, options);


function theta_min = find_theta_star(Q, W, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    coarse_res = 100;
    theta = linspace(min(W), max(W), coarse_res);
    derivative = tp_derivative(theta, Q, tt, hat_phi_W, sqrt_psi_hat_W, weight);
    [~, I] = min(derivative);
    
    I_low = max(I - 1, 1);
    I_high = min(I + 1, coarse_res);
    
    %     figure(1)
    %     plot(theta, derivative)
    %     drawnow;
    
    fine_res = 100;
    theta = linspace(theta(I_low), theta(I_high), fine_res);
    derivative = tp_derivative(theta, Q, tt, hat_phi_W, sqrt_psi_hat_W, weight);
    [~, I] = min(derivative);
    theta_min = theta(I);
end

function Q = exchange_mass(Q, theta_star, I, objective_func)
    coarse_res = 100;
    p_coarse = linspace(0, 1, coarse_res);
    tp_coarse = linspace(1, coarse_res);
    
    for i = 1:coarse_res
        p_test = p_coarse(i);
        pj = Q.ProbWeights;
        pj = shift_mass(pj, p_test, theta_star, I);
        
        % Calculate tp at new distribution
        tp_coarse(i) = objective_func(pj);
    end
    
    [~, I_min] = min(tp_coarse);
    I_low = max(I_min - 1, 1);
    I_high = min(I_min + 1, coarse_res);
    fine_res = 100;
    p_fine = linspace(p_coarse(I_low), p_coarse(I_high), fine_res);
    tp_fine = linspace(1, fine_res);
    
    for i = 1:fine_res
        p_test = p_fine(i);
        pj = Q.ProbWeights;
        pj = shift_mass(pj, p_test, theta_star, I);
        
        % Calculate tp at new distribution
        tp_fine(i) = objective_func(pj);
    end
    
    [~, I_star] = min(tp_fine);
    p_star = p_fine(I_star);
    
%     plot(p_coarse, tp_coarse);
%     p_star
    
    Q.ProbWeights = shift_mass(Q.ProbWeights, p_star, theta_star, I);
    
    remove_mass = Q.ProbWeights <= 0.000001;
    Q.Support(remove_mass) = [];
    Q.ProbWeights(remove_mass) = [];
    Q.ProbWeights = Q.ProbWeights / sum(Q.ProbWeights);
end

function pj = shift_mass(pj, p_star, theta_star, I)
    pj(theta_star) = p_star; %Put mass p_test at theta_star
        
    %Remove mass p_test from the prob weights in reverse order from I
    %(ie move mass from worst masses to the best mass)
    freddy = (cumsum(pj(I(end:-1:1)))) < p_star;
    pj(I(freddy(end:-1:1))) = 0;

    % Partial amount to remove from 
    leftover = sum(pj) - 1;
    % Indices we haven't removed stuff from 
    part_index = I(~freddy(end:-1:1));
    % Last index is the one we want to remove some things from
    part_index = part_index(end);
    pj(part_index) = pj(part_index) - leftover;
end