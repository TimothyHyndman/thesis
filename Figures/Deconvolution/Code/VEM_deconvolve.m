function [Q, tt, optim_values_VEM] = VEM_deconvolve(W, xx, Q_init)
    
    % Initial distribution
    if nargin < 3
        Q.Support = mean(W);
        Q.ProbWeights = 1;
    else
        Q = Q_init;
    end

    % VEM like algorithm for TP
    use_penalties = true;
    n = length(W);
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

    objective_func = @(pj, xj) tp_and_penalties(tt, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight, use_penalties);
    % -------------------------------------------------------------------------

    optim_values_VEM.tp_initial = calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight);
    [penalty1, penalty2] = penalties(Q.ProbWeights, Q.Support, tt, hat_phi_W);
    optim_values_VEM.penalty1_initial = penalty1;
    optim_values_VEM.penalty2_initial = penalty2;
    optim_values_VEM.objective_func_initial = objective_func(Q.ProbWeights, Q.Support);
    optim_values_VEM.var_initial = calculate_var(Q);

    looping = true;
    counter = 0;
    tp_old = Inf;
    tp_update_tol = 0;
    max_n_loops = 1000;

    while looping
        % Add theta which minimizes tp to support
        theta_min = find_theta_star(Q, W, objective_func);
        Q.Support = [Q.Support, theta_min];
        Q.ProbWeights = [Q.ProbWeights, 0];

        %Rank support points by their tp_derivative value
        theta_derivatives = tp_derivative(Q.Support, Q, objective_func);
        [~, I] = sort(theta_derivatives);

        %Move mass from highest tp_masses to lowest
    %     objective_func = @(pj) calculate_tp(tt, pj, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight);
        objective_func_2 = @(pj) tp_and_penalties(tt, pj, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight, use_penalties);
        Q = exchange_mass(Q, I(1), I, objective_func_2);

        % Merge points that are close together
        Q = merge_points(Q, objective_func);

        % Check for convergence
        tp_new = objective_func(Q.ProbWeights, Q.Support);
        if abs(tp_new - tp_old) <= tp_update_tol
            looping = false;
        else
            tp_old = tp_new;
        end

        counter = counter + 1;
        if counter > max_n_loops
            looping = false;
        end
    end

    optim_values_VEM.n_loops = counter - 1;
    optim_values_VEM.tp_final = calculate_tp(tt, Q.ProbWeights, Q.Support, hat_phi_W, sqrt_psi_hat_W, weight);
    [penalty1, penalty2] = penalties(Q.ProbWeights, Q.Support, tt, hat_phi_W);
    optim_values_VEM.penalty1_final = penalty1;
    optim_values_VEM.penalty2_final = penalty2;
    optim_values_VEM.objective_func_final = objective_func(Q.ProbWeights, Q.Support);
    optim_values_VEM.var_final = calculate_var(Q);
    
    optim_values_VEM = orderfields(optim_values_VEM);

end

function yy = var_derivative(theta, Q)
    xj = Q.Support;
    pj = Q.ProbWeights;
    mu = pj*xj';
    
    yy = (theta - mu).^2 - pj*(xj - mu).^2.';
end

function sigma = calculate_var(Q)
    mu = Q.ProbWeights(:)' * Q.Support(:);
    
    sigma = Q.ProbWeights(:)' * (Q.Support(:) - mu).^2;
end

function Q = merge_points(Q, objective_func)
    [xj, I] = sort(Q.Support);
    pj = Q.ProbWeights(I);
    if length(xj) == 1
        looping = false;
    else
        looping = true;
    end
    index = 1;
    while looping
        
        tp_current = objective_func(pj, xj);
        
        % Try combining points index and index+1
        xj_test = xj;
        pj_test = pj;
        xj_test(index+1) = (pj(index) * xj(index)+pj(index+1)*xj(index+1))/(pj(index)+pj(index+1));
        pj_test(index+1) = pj(index)+pj(index+1);
        pj_test(index) = [];
        xj_test(index) = [];
        
        tp_test = objective_func(pj_test, xj_test);
        
        if tp_test <= tp_current
            xj = xj_test;
            pj = pj_test;
        else
            index = index + 1;
        end
        
        if index >= (length(xj) - 1)
            looping = false;
        end
        
    end
    Q.Support = xj;
    Q.ProbWeights = pj;
    
end

function theta_min = find_theta_star(Q, W, objective_func)
    coarse_res = 100;
    theta = linspace(min(W), max(W), coarse_res);
    derivative = tp_derivative(theta, Q, objective_func);
    [~, I] = min(derivative);
    
    I_low = max(I - 1, 1);
    I_high = min(I + 1, coarse_res);
    
    figure(1)
    plot(theta, derivative)
    drawnow;
    
    fine_res = 20;
    for j = 1:5
        theta = linspace(theta(I_low), theta(I_high), fine_res);
        derivative = tp_derivative(theta, Q, objective_func);
        [~, I] = min(derivative);
        I_low = max(I - 1, 1);
        I_high = min(I + 1, fine_res);
    end
    
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
    
    
    fine_res = 20;
    for j = 1:5
        p_fine = linspace(p_coarse(I_low), p_coarse(I_high), fine_res);
        tp_fine = zeros(1, fine_res);

        for i = 1:fine_res
            p_test = p_fine(i);
            pj = Q.ProbWeights;
            pj = shift_mass(pj, p_test, theta_star, I);

            % Calculate tp at new distribution
            tp_fine(i) = objective_func(pj);
        end
        
        [~, I_min] = min(tp_fine);
        I_low = max(I_min - 1, 1);
        I_high = min(I_min + 1, fine_res);
        
        p_coarse = p_fine;
    end
    
    p_star = p_fine(I_min);
    
%     plot(p_coarse, tp_coarse);
%     p_star
    
    Q.ProbWeights = shift_mass(Q.ProbWeights, p_star, theta_star, I);
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

function [penalty1, penalty2] = penalties(pj, xj, t, hat_phi_W)

    re_hat_phi_W = real(hat_phi_W);
    im_hat_phi_W = imag(hat_phi_W);
    norm_hat_phi_W = sqrt(re_hat_phi_W.^2 + im_hat_phi_W.^2);
    [re_phi_p, im_phi_p, norm_phi_p] = computephiX(t, xj, pj);

    %Need phi_U to be real
    a = re_phi_p(:) .* im_hat_phi_W(:);
    b = im_phi_p(:) .* re_hat_phi_W(:);
    penalty1  = sum(abs(a - b));

    %impose a penalty if |phi_U| is greater than 1:
    hat_phi_U = norm_hat_phi_W(:) ./ norm_phi_p(:);
    penalty2 = sum(hat_phi_U(hat_phi_U > 1));
end

function fval = tp_and_penalties(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight, use_penalties)

    fval = calculate_tp(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight);
    
    if use_penalties
        [penalty1, penalty2] = penalties(pj, xj, t, hat_phi_W);
        fval = fval + 500*penalty1 + 500*penalty2;
    end
end