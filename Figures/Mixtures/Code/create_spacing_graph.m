%spacing graph

res = 50;
spacing = linspace(0,2,res);

sigma = 1;
phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));

n = 10;
m = zeros(1,res);

for i = 1:res
	a = spacing(i);
	X = linspace(0, (n-1)*a, n)';

	Q = VEM(phi,X);

	m(i) = number_points_support(phi, X, Q);
    figure(1)
    plot(spacing, m)
    drawnow
end

plot(spacing, m)