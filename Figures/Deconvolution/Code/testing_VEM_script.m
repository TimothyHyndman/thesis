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

% U
sigU = sqrt(NSR * varX);
U = normrnd(0, sigU, 1, n);

W = X + U;
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
    Q = exchange_mass(Q, I, objective_func);
    
    counter = counter + 1;
    if counter > 100
        looping = false;
    end
end

calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight)

figure(2)
scatter(Q.Support, Q.ProbWeights)
drawnow


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

function Q = exchange_mass(Q, I, objective_func)
    coarse_res = 100;
    p_coarse = linspace(0, 1, coarse_res);
    tp_coarse = linspace(1, coarse_res);
    
    for i = 1:coarse_res
        p_test = p_coarse(i);
        pj = Q.ProbWeights;
        pj = shift_mass(pj, p_test, I);
        
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
        pj = shift_mass(pj, p_test, I);
        
        % Calculate tp at new distribution
        tp_fine(i) = objective_func(pj);
    end
    
    [~, I_star] = min(tp_fine);
    p_star = p_fine(I_star);
    
%     plot(p_coarse, tp_coarse);
%     p_star
    
    Q.ProbWeights = shift_mass(Q.ProbWeights, p_star, I);
end

function pj = shift_mass(pj, p_star, I)
    pj(I(1)) = p_star; %Put mass p_test at theta_min
        
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