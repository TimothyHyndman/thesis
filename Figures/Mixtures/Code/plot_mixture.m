function fig = plot_mixture(phi, Q, Y, options, Q_unsimplified)

    options = fill_default_options(options);

    %Colours
    % highlight = [255 189 76]/255;
    highlight = [0 51 104]/255;
    primary = [38 82 127]/255;
    primarylight = [102,177,255]/255;
    white = [255,255,255]/255;
    backgroundblue = 0.8*white+0.2*primarylight;
    grey = [150, 150, 150]/255;


    % Mixture Density
    xx = linspace(min(Y)-3,max(Y)+3,1000);
    yy = 0*xx;
    for i = 1:length(Q.ProbWeights)
        yy = yy + Q.ProbWeights(i)*phi(xx',Q.Support(i))';
    end

    fig = figure('pos',[200 200 800 450]);
    legend();

    if options.plot_hist
        hold on
        h_histogram = histogram(Y,'Normalization','pdf','FaceColor',primary,'FaceAlpha',1);
    end

    if options.plot_data_points
        hold on
        h_datapoints = scatter(Y,0*Y,'*')
    end

    if options.plot_mixture_density
        hold on
        h_density = plot(xx,yy,'m','Color',highlight,'LineWidth',4, 'DisplayName', 'f_Q(x)');
    end

    if options.plot_Q_unsimplified
        hold on
        h_Qunsimplified = scatter(Q_unsimplified.Support,Q_unsimplified.ProbWeights,'filled','MarkerFaceColor',grey, 'DisplayName', 'Q*');
    end

    if options.plot_Q
        hold on
        h_Q = scatter(Q.Support,Q.ProbWeights,'filled','MarkerFaceColor',highlight,'DisplayName', 'Q');
    end
   
    hold off

    %Background Colour
    % set(gca,'Color',backgroundblue);

    xlabel('x','FontSize',20)

    fig = gca;

    if options.save
    %Print to file... Make sure that the current folder is where you want to
    %save it!!!
    if strcmp(options.filename, '')
        filename = 'no_filename_set_mixture_plot.png';
    else
        filename = options.filename;
    end

    saveas(fig, filename)
    end
end

function options = fill_default_options(options)
    if ~isfield(options, 'plot_hist')
        options.plot_hist = false;
    end

    if ~isfield(options, 'plot_data_points')
        options.plot_data_points = false;
    end

    if ~isfield(options, 'plot_mixture_density')
        options.plot_mixture_density = true;
    end

    if ~isfield(options, 'plot_Q')
        options.plot_Q = false;
    end

    if ~isfield(options, 'plot_Q_unsimplified')
        options.plot_Q_unsimplified = false;
    end

    if ~isfield(options, 'save')
        options.save = false;
    end

    if ~isfield(options, 'filename')
        options.filename = '';
    end
end