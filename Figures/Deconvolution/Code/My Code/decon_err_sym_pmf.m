%% Deconvolution of data X

function [Q, tt, normhatphiW, optim_values] = decon_err_sym_pmf(W, m, n_tp_iter, n_var_iter, show_diagnostics, decon_options)
    n = length(W);
    
    decon_options = fill_default_options(decon_options);
    

    diagnostic = @(message) print_diagnostic(message, show_diagnostics);
    
    % Precalculate phi_W -------------------------------------------------------
    length_tt = 100;
    mu_K2 = 6;
    RK = 1024 / 3003 / pi;
    hnaive = ((8 * sqrt(pi) * RK/3/mu_K2^2)^0.2) * sqrt(var(W)) * n^(-1/5);
    hmin = hnaive/3;
    tt = linspace(-1/hmin, 1/hmin, length_tt);

    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    normhatphiW = sqrt(rehatphiW.^2 + imhatphiW.^2)';
    t_star = find_t_cutoff(normhatphiW, tt, n);
    
%     plot(tt, normhatphiW)
%     t_star
    
    
%     tmp=tt(normhatphiW<n^(-0.25));
%     tt1=max(tmp(tmp<0))
    
    
    tt = linspace(-t_star, t_star, length_tt);

    [rehatphiW, imhatphiW] = compute_phi_W(tt, W);
    hat_phi_W = complex(rehatphiW, imhatphiW).';
    [~, ~, sqrt_psi_hat_W] = compute_psi_W(tt, W);
    sqrt_psi_hat_W = sqrt_psi_hat_W';

    weight = KernelWeight('Epanechnikov',tt);

    %--------------------------------------------------------------------------%
    % Solve optimization problem to find PMF
    %--------------------------------------------------------------------------%

    options = optimoptions('fmincon', ...
                           'Display', 'off', ...
                           'Algorithm', 'active-set', ...
                           'TolFun', 1e-6, ...
                           'MaxFunctionEvaluations', 1e4, ...
                           'MaxIterations', 5e3);
    
    [A, B] = create_bound_matrices(W, m);

    % ------------------
    % Min T(p)
    % ------------------
    diagnostic("Minimizing T(p)")

    func = @(x) tp_objective(x, m, tt, hat_phi_W, sqrt_psi_hat_W, weight, decon_options);

    fmax = Inf;
    counter = 0;
    while counter < n_tp_iter
        pj_0 = unifrnd(0, 1, [1, m]);
        pj_0 = pj_0 / sum(pj_0);
        xj_0 = sort(unifrnd(min(W), max(W), [1, m]));

        x0 = [pj_0(1:end-1), xj_0]';
        [x, fval, exitflag] = fmincon(func, x0, A, B,[],[],[],[],[], options);
        [pj_new, xj_new] = x_to_pmf(x);

        diagnostic(num2str(exitflag))

        if fval < fmax && exitflag >= 0
            fmax = fval;
            pj = pj_new;
            xj = xj_new;
            diagnostic(num2str(fval))
        end

        if (exitflag >= 0)
            counter = counter + 1;
        end
    end
    
    % USE VEM TO GET BETTER RESULT
%     Q_tp.Support = xj;
%     Q_tp.ProbWeights = pj;
%     [Q, ~, optim_values_VEM] = VEM_deconvolve(W, Q_tp);
%     
%     optim_values_VEM
%     pj = Q.ProbWeights;
%     xj = Q.Support;

    % ------------------
    % Min Var
    % ------------------

    %Calculate penalties once 
    tp_max = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    [penalty1_max, penalty2_max, ~] = penalties(pj, xj, tt, hat_phi_W);
    diagnostic(join(["tp_max =", num2str(tp_max)]))
    diagnostic(join(["penalties =", num2str(penalty1_max), ",", num2str(penalty2_max)]))

    %rescale penalties to allow some room
    penalty_tolerance_scale = 0.05;
    penalty_tolerance_scale = 0;
    tp_max = tp_max * (1 + penalty_tolerance_scale);
    penalty1_max = penalty1_max * (1 + penalty_tolerance_scale);
    penalty2_max = penalty2_max * (1 + penalty_tolerance_scale);

    %Find initial value for varmin based on best solution for fmax
    varmin = var_objective([pj(1:end-1), xj]');
    varmin_init = varmin;

    func = @(x) var_objective(x);
    nonlcon = @(x)phaseconstraint(x, tp_max,penalty1_max, penalty2_max,tt,hat_phi_W,sqrt_psi_hat_W,weight, decon_options);

    counter = 0;
    
    diagnostic("Minimizing Variance")
    first_try = true;
    while counter < n_var_iter
        if first_try
            pj_0 = pj;
            xj_0 = xj;
            first_try = false;
        else
            pj_0 = unifrnd(0,1,[1,m]);
            pj_0 = pj_0 / sum(pj_0);
            xj_0 = sort(unifrnd(min(W), max(W), [1,m]));
        end

        x0 = [pj_0(1:end-1), xj_0]';
        [x, fval, exitflag] = fmincon(func, x0, A, B, [], [], [], [], nonlcon, options);
        [pj_new, xj_new] = x_to_pmf(x);

        if ~is_feasible(pj_new, xj_new, A, B, nonlcon)
            exitflag = -99;
%             exitflag = 1;
        end

        diagnostic(num2str(exitflag))

        if fval < varmin && exitflag >= 0
            varmin = fval;
            pj = pj_new;
            xj = xj_new;
            diagnostic(num2str(fval))
        end

        if (exitflag >= 0)
            counter = counter + 1;
        else
            counter = counter;
        end

        drawnow
    end

    tp_final = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    [penalty1_final, penalty2_final, ~] = penalties(pj, xj, tt, hat_phi_W);
    var_final = var_objective([pj(1:end-1), xj]');

    diagnostic(join(["Initial variance was", num2str(varmin_init)]))
    diagnostic(join(["Final variance is", num2str(var_final)]))
    diagnostic(join(["T(p) =", num2str(tp_final)]))
    diagnostic(join(["penalties =", num2str(penalty1_final), ",", num2str(penalty2_final)]))

    % Finalize -----------------------------------------------------------------
%     [pj,xj] = simplify_masses(pj,xj);
    Q.Support = xj;
    Q.ProbWeights = pj;
    
    optim_values.tp_final = tp_final;
    optim_values.var_final = var_final;
    optim_values.penalty1_final = penalty1_final;
    optim_values.penalty2_final = penalty2_final;
    optim_values.objective1_final = tp_final + 500*(penalty1_final + penalty2_final);
end

%-------------------------------------------------------------------------------
% Local Functions
%-------------------------------------------------------------------------------

function [pj, xj] = x_to_pmf(x)
    m = (length(x) + 1) / 2;
    x = x';
    pj = [x(1:m-1), 1 - sum( x(1:m-1))];
    xj = x(m:(2 * m - 1));
end

function flag = is_feasible(pj, xj, A, B, nonlcon)
    flag = 1;
    x = [pj(1:end-1), xj]';

    lin_tol = 1e-5;
    if any(A * x - B > lin_tol)
        flag = 0;
    end

    non_lin_tol = 1e-5;
    if any(nonlcon(x) > non_lin_tol)
        flag = 0;
    end
end

function [fval, penalty1, penalty2, tp] = tp_objective(x, m, tt, hat_phi_W, sqrt_psi_hat_W, weight, decon_options)
    %Extract weights and roots
    [pj, xj] = x_to_pmf(x);

    %Integrate
    tp = calculate_tp(tt,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);

    %Add penalty terms
    if decon_options.penalties
        [penalty1, penalty2, ~] = penalties(pj,xj,tt,hat_phi_W);
    else
        penalty1 = 0;
        penalty2 = 0;
    end
    
    penalty_scale = 500;
    fval = tp + penalty_scale * (penalty1 + penalty2);
end

function var1 = var_objective(x)
    [pj, xj] = x_to_pmf(x);
    mean1 = sum(pj * xj');
    var1 = sum(pj * (xj - mean1)'.^2);
end

function [c,ceq] = phaseconstraint(x, tp_max, penalty1_max, penalty2_max, t, hat_phi_W, sqrt_psi_hat_W, weight, decon_options)
    [pj, xj] = x_to_pmf(x);
    tp = calculate_tp(t,pj,xj,hat_phi_W,sqrt_psi_hat_W,weight);
    
    if decon_options.penalties
        [penalty1, penalty2, ~] = penalties(pj,xj,t,hat_phi_W);
    else
        penalty1 = 0;
        penalty2 = 0;
    end

    c = [tp - tp_max, penalty1 - penalty1_max, penalty2 - penalty2_max];
%     c = [tp - tp_max, penalty1 - penalty1_max];    %This line matches Aurore's original code (except she just has penalty1 - 0)
    ceq=[];
end

function [penalty1, penalty2, penalty3] = penalties(pj, xj, t, hat_phi_W)
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
    penalty2 = sum(hat_phi_U(hat_phi_U > 1) - 1);
    
    penalty3 = 0;   %Removed this penalty
end

function [pj, xj] = simplify_masses(pj, xj)
    zero_tolerance = 0.001;
    search_size = 0.001;

    %Get rid of zeros
    index = 1;
    looping = true;
    while looping
        if pj(index) < zero_tolerance
            xj(index:end-1) = xj(index+1:end);
            xj(end) = [];
            pj(index:end-1) = pj(index+1:end);
            pj(end) = [];
        else
            index = index+1;
        end
        if index >= length(xj)
            looping = false;
        end
    end

    if pj(end) < zero_tolerance
        pj(end) = [];
        xj(end) = [];
    end

    %Get rid of duplicates
    index = 1;
    looping = true;
    if length(xj) < 2
        looping = false;
    end
    while looping
        if abs(xj(index) - xj(index+1)) < search_size
            xj(index) = (xj(index)*pj(index)+xj(index+1)*pj(index+1))/(pj(index)+pj(index+1));  %Weighted average location
            xj(index+1:end-1) = xj(index+2:end);
            pj(index) = pj(index)+pj(index+1);
            pj(index+1:end-1) = pj(index+2:end);
            pj(end) = [];
            xj(end) = [];
        else
            index = index+1;
        end
        if index == length(xj)
            looping = false;
        end
    end

    pj = pj/sum(pj);    %Normalise
end

function [A, B] = create_bound_matrices(W, m)
    % pj non-negative
    A = zeros(2*m + 1, 2*m - 1);
    A(1:(m-1), 1:(m-1)) = eye( m - 1 );
    B = zeros(2 * m - 1, 1);
    % pj sum to less than 1
    A(m, 1:m-1) = -ones(1,m-1);
    B(m) = -1;

    % thetaj are increasing
    for i = 1:(m-1)
        A(m+i, (m+i-1):(m+i)) = [-1, 1];
    end
    % min(W) < thetaj < max(W)
    A(2*m, m) = 1;
    A(2*m+1, 2*m-1) = -1;
    B(2*m) = min(W);
    B(2*m+1) = -max(W);

    A = -A;
    B = -B;
end

function print_diagnostic(message, show_diagnostics)
    if show_diagnostics
        disp(message)
        drawnow
    end
end

function options = fill_default_options(options)
    if ~isfield(options, 'penalties')
        options.penalties = true;
    end
end