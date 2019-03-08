%Plot number of components required as phase diagram
res = 512;   %ultra-1024, vhigh-512, high-256, temp-32
width = 5;


%Normal
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));

Z_norm = calculate_flag_graph(phi, width, res);

options.line_width = 1;
options.filename = 'normal_flag_graph.png';
plot_flag_graph(Z_norm, sigma, width, options);


%Plot number of components required as phase diagram

%Cauchy
lambda = sqrt(3);   %makes inflection points at +- 1
phi = @(x,u) (1/(pi*lambda))*(1 + ((repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))./lambda).^2).^(-1);

Z_cauchy = calculate_flag_graph(phi, width, res);

options.line_width = 1;
options.filename = 'cauchy_flag_graph.png';
plot_flag_graph(Z_cauchy, sigma, width, options);


figure()

options.cmap = 'parula';
options.cmap = [0 51 104; 19 67 116; 38 82 127; 57 123 191; 76 164 255]/255;
options.filename = 'norm_cauchy_compare_flag_graph.png';
options.save = false;
plot_flag_graph((Z_norm+Z_cauchy)/2, sigma, width, options)