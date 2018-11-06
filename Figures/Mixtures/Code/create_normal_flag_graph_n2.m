%Plot number of components required as phase diagram

%Normal
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));


res = 512;   %ultra-1024, vhigh-512, high-256, temp-32
width = 5;

Z = calculate_flag_graph_n2(phi, width, res);


options.line_width = 1;
options.save = true;
options.filename = 'normal_flag_graph_n2.png';
options.n = 2;
plot_flag_graph(Z, sigma, width, options);