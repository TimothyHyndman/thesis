% T(p) derivative


function derivative = tp_derivative(theta, Q, objective_func)   
    xj = Q.Support;
    pj = Q.ProbWeights;
    
    delta_pj = 1e-10;
    derivative = zeros(1, length(theta));
    
%     current_tp = calculate_tp(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight);    
    current_tp = objective_func(pj, xj);
    
    for i = 1:length(theta)
        xj_new = [xj, theta(i)];
        pj_new = [pj, delta_pj];
        pj_new=  pj_new / sum(pj_new);
        
%         derivative(i) = (calculate_tp(t, pj_new, xj_new, hat_phi_W, sqrt_psi_hat_W, weight) - current_tp)/delta_pj; 
        derivative(i) = (objective_func(pj_new, xj_new) - current_tp)/delta_pj; 
    end
end