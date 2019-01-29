% Save deconvolution graph
function fig = plot_deconvolution_graph(Q, xx, yy, W, options)

    options = fill_default_options(options);
    
    highlight = [255 189 76]/255;
    primarydark = [0 51 104]/255;
    primary = [38 82 127]/255;
    primarylight = [102,177,255]/255;
    white = [255,255,255]/255;
    backgroundblue = 0.8*white+0.2*primarylight;
    grey = [150, 150, 150]/255;
    
    fig = figure('pos',[200 200 800 450]);
    
    hold on
    
    if options.plot_histogram
        fig_histogram = histogram(W,'Normalization','pdf', 'FaceColor', primary,'FaceAlpha',1, 'DisplayName','W');
    end
    
    if options.plot_true_dens
        fig_true_dens = plot(xx, options.true_dens(xx), 'Color', primarylight,'LineWidth', options.line_width, 'DisplayName', 'True f_X');
    end
    
    if options.plot_true_pmf
        fig_true_pmf = scatter(options.true_pmf.Support, options.true_pmf.ProbWeights, 'filled','MarkerFaceColor', primarylight, 'DisplayName', 'True f_X');
    end
    
    if options.plot_naive
        n = length(W);
        h = options.naivebw;
        xout=outerop(xx,W,'-');
        fWEF=normpdf(xout,0,h)*ones(n,1)/n;
        fig_naive = plot(xx, fWEF, ':', 'Color', primarylight,'LineWidth', options.line_width, 'DisplayName', 'Naive f_X');
    end

    if options.plot_density
        fig_density = plot(xx, yy, 'm','Color', primarydark,'LineWidth', options.line_width, 'DisplayName', 'Estimated f_X');
    end
    
    if options.plot_masses
        fig_masses = scatter(Q.Support, Q.ProbWeights, 'filled','MarkerFaceColor', grey, 'DisplayName', 'Estimated F_X');
    end
    
    hold off
    
    xlabel('x', 'FontSize', 20)
    
    yl = ylim;
    ylim([0, yl(2)])
    
    legend('show')
    
%     fig = gca;
    
    if options.save
        %Print to file... Make sure that the current folder is where you want to
        %save it!!!
        if strcmp(options.filename, '')
            filename = 'no_filename_set_deconvolution_plot.png';
        else
            filename = options.filename;
        end
        
        saveas(fig, filename)
    end
end

function options = fill_default_options(options)
    if ~isfield(options, 'save')
        options.save = false;
    end
    if ~isfield(options, 'line_width')
        options.line_width = 4;
    end
    if ~isfield(options, 'filename')
        options.filename = '';
    end
    if ~isfield(options, 'plot_histogram')
        options.plot_histogram = false;
    end
    if ~isfield(options, 'plot_density')
        options.plot_density = false;
    end
    if ~isfield(options, 'plot_masses')
        options.plot_masses = true;
    end
    if ~isfield(options, 'true_dens')
        options.true_dens = @(x) 0;
    end
    if ~isfield(options, 'plot_true_dens')
        options.plot_true_dens = false;
    end
    if ~isfield(options, 'true_pmf')
        options.true_pmf = [];
    end
    if ~isfield(options, 'plot_true_pmf')
        options.plot_true_pmf = false;
    end
    if ~isfield(options, 'plot_naive')
        options.plot_naive = false;
    end
    if ~isfield(options, 'naivebw')
        options.naivebw = NA;
    end
end