function [yy, Q, tt, optim_values] = moving_mass_deconvolve(W, xx)
    
   [yy, Q, tt, optim_values] = decon_err_sym(W, xx, 10);

end