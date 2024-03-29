function tp = calculate_tp(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight)
%     t = t(:);
    hat_phi_W = hat_phi_W(:);
    sqrt_psi_hat_W = sqrt_psi_hat_W(:);
    weight = weight(:);
    
    %Calculate characteristic function of our discrete distribution
    pj = pj(:).';
    xj = xj(:);
    t = t(:);
    [re_phi_p, im_phi_p, norm_phi_p] = computephiX(t, xj, pj);
    phi_p = complex(re_phi_p, im_phi_p);
    
    norm_phi_p = norm_phi_p(:);
    phi_p = phi_p(:);

    %Calculate integrand    
    % This is what Aurore uses which matches the version in the paper if
    % you read it properly
    integrand = abs(hat_phi_W .* norm_phi_p - sqrt_psi_hat_W .* phi_p).^2.*weight;
    dt = t(2) - t(1);
    tp = dt * sum(integrand);
end