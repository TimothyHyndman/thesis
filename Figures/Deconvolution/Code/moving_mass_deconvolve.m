function [yy, Q, tt, optim_values] = moving_mass_deconvolve(W, xx, m, decon_options)
   
tic
    if nargin < 3
        m = 10;
    end
   [yy, Q, tt, optim_values] = decon_err_sym(W, xx, m, decon_options);
   optim_values = orderfields(optim_values);
toc
end