% T(p) derivative


function tp_derivative(W, Q, t, hat_phi_W, sqrt_psi_hat_W, weight)
    res = 100;
    theta = linspace(min(W), max(W), res);
    
    xj = Q.Support;
    pj = Q.ProbWeights;
    
    delta_pj = 1e-12;
    derivative = zeros(1, res);
    
    current_tp = calculate_tp1(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight);
    
    for i = 1:res
        xj_new = [xj, theta(i)];
        pj_new = [pj, delta_pj];
        pj_new=  pj_new / sum(pj_new);
        
        derivative(i) = (calculate_tp1(t, pj_new, xj_new, hat_phi_W, sqrt_psi_hat_W, weight) - current_tp)/delta_pj; 
    end
    
    
    plot(theta, derivative)
end