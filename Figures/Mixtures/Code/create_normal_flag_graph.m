%Plot number of components required as phase diagram

%Normal
sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));
%Uniform
% sigma = 1;
% phi = @(x,u) (1/2*sigma)*(abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))<=sigma);
%Triangle
% sigma = 1;
% phi = @(x,u) (1/sigma)*(1 - (1/sigma)*abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))).*(abs(repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))<=sigma);
%Cauchy
% lambda = 1;
% phi = @(x,u) (1/(pi*lambda))*(1 + ((repmat(x,[1,length(u)]) - repmat(u,[length(x),1]))./lambda).^2).^(-1);


res = 1024;   %ultra-1024, vhigh-512, high-256
width = 5;

Z = calculate_flag_graph(phi, width, res);


options.line_width = 1;
options.save = true;
options.filename = '';
plot_flag_graph(Z, sigma, width, options);