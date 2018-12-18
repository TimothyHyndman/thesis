function tp = calculate_tp1(t, pj, xj, hat_phi_W, sqrt_psi_hat_W, weight)

    %Calculate characteristic function of our discrete distribution
    [re_phi_p, im_phi_p, norm_phi_p] = computephiX(t, xj, pj);
    phi_p = complex(re_phi_p, im_phi_p);

    %Calculate integrand
    integrand = abs(hat_phi_W - sqrt_psi_hat_W .* phi_p ./ norm_phi_p).^2.*weight';
    dt = t(2) - t(1);
    tp = dt * sum(integrand);
end