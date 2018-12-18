function t_star = find_t_cutoff(normhatphiW, tt, n)
    ind = find(tt >= 0);
    d_phi_W = normhatphiW(ind(2:end)) - normhatphiW(ind(1:end-1));

    if isempty(find(d_phi_W >= 0))
        t_star = tt(end);
    else
        first_min_ind = ind(min(find(d_phi_W >= 0)));
        phi_W_threshold = max(normhatphiW(ind(ind >= first_min_ind)));
        tmp = tt(normhatphiW <= phi_W_threshold);
        t_star = min(tmp(tmp>0));
    end
    
    % Alternative method (matches Aurore's code)        
    tmp=tt(normhatphiW<n^(-0.25));
    tt1=max(tmp(tmp<0));
    t_star = abs(tt1);
end