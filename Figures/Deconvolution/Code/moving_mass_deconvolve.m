function [yy, Q, tt, optim_values] = moving_mass_deconvolve(W, xx, m)
   
    if nargin < 3
        m = 10;
    end
   [yy, Q, tt, optim_values] = decon_err_sym(W, xx, m);
   optim_values = orderfields(optim_values);

end