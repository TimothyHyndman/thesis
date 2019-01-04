function [fXdeconvoluted, Q, tt, optim_values] = decon_err_sym(W, xx, m, pmf, bw, show_diagnostics)

    if (~exist('xx','var') | isempty(xx))
        xx = linspace(min(W), max(W), 100);
    end

    if (~exist('m', 'var') | isempty(m))
      m = 10;
    end

    if m < 2
      error("m must be at least 2")
    end

    if (~exist('show_diagnostics', 'var') | isempty(show_diagnostics))
      show_diagnostics = 0;
    end

    if (~exist('pmf', 'var') | isempty(pmf))
      pmf = 0;
    end

    if (~exist('bw', 'var') | isempty(bw))
      bw = [];
    end


    % Deconvolve to pmf --------------------------------------------------------
    n_tp_iter = 20;
    n_var_iter = 20;
    [Q, tt, normhatphiW, optim_values] = decon_err_sym_pmf(W, ...
                                             m, ...
                                             n_tp_iter, ...
                                             n_var_iter,...
                                             show_diagnostics);

    % Convert pmf to pdf -------------------------------------------------------   
    if pmf
      fXdeconvoluted = [];
    else
      fXdeconvoluted = decon_err_sym_pmf2pdf(xx, ...
                                             tt, ...
                                             Q.Support, ...
                                             Q.ProbWeights, ...
                                             W, ...
                                             bw);
    end
end