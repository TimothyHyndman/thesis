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

looping = true;
counter = 0;
while looping
    [theta_min, theta_max] = find_theta_star(Q, W, tt, hat_phi_W, sqrt_psi_hat_W, weight);
    p_min = find_p_min(Q, theta_min, tt, hat_phi_W, sqrt_psi_hat_W, weight);

    Q.Support = [Q.Support, theta_min];
    Q.ProbWeights = [Q.ProbWeights * (1 - p_min), p_min]; 
    I_remove = Q.ProbWeights <= 0;
    Q.Support(I_remove) = [];
    Q.ProbWeights(I_remove) = [];
    
    counter = counter + 1;
    if counter > 100
        looping = false;
    end
end


calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight)

function [theta_min, theta_max] = find_theta_star(Q, W, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    res = 1000;
    theta = linspace(min(W), max(W), res);

    derivative = tp_derivative(theta, Q, tt, hat_phi_W, sqrt_psi_hat_W, weight);

    figure(1)
    plot(theta, derivative)
    drawnow;

    [~, I] = min(derivative);
    theta_min = theta(I);
    
    [~, I] = max(derivative);
    theta_max = theta(I);
end

function p_min = find_p_min(Q, theta_star, tt, hat_phi_W, sqrt_psi_hat_W, weight)
    p_res = 1000;
    p_test = linspace(0, 1, p_res);
    tp_values = zeros(1, p_res);

    for i = 1:p_res
        xj_test = [Q.Support, theta_star];
        pj_test = [Q.ProbWeights * (1 - p_test(i)), p_test(i)];
        tp_values(i) = calculate_tp(tt, pj_test, xj_test, hat_phi_W, sqrt_psi_hat_W, weight);
    end

    [~, I] = min(tp_values);
    p_min = p_test(I);
end