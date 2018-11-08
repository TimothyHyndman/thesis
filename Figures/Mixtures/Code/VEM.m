% X is column vector

function [Q,log_likelihood] = VEM(phi,X)

	%Initial value
	Q.Support = mean(X);
	Q.ProbWeights = 1;

	looping = true;

	while looping
		% Find theta* = argmax D(theta)
		[theta_star, derivative] = max_theta2(phi, X, Q);

		if derivative < 1e-3
			looping = false;
		end

		% Add theta* to support
		Q.Support = [Q.Support, theta_star];
		Q.ProbWeights = [Q.ProbWeights, 0];

		% Find point of support theta- which minimizes D(theta_j)
		[~, index_min_support] = min(D(phi, X, Q, Q.Support));

		% Find how much mass to move from theta- to theta*
		fun = @(x) -blah(phi, X, Q, index_min_support, x);

		res = 100;
		xx = linspace(0, Q.ProbWeights(index_min_support), res);
		for i = 1:res
			likelihoods(i) = -fun(xx(i));
		end
		figure(1)
		plot(xx, likelihoods)
		drawnow
		delta = fminbnd(fun, 0, Q.ProbWeights(index_min_support))

		% res = 100;
		% h = linspace(0, Q.ProbWeights(index_min_support), res);
		% likelihoods = zeros(1, res);
		% for i = 1:res
		% 	Q_test = Q;
		% 	Q_test.ProbWeights(end) = h(i);
		% 	Q_test.ProbWeights(index_min_support) = Q_test.ProbWeights(index_min_support) - h(i);
		% 	likelihoods(i) = L(phi, X, Q_test);
		% end
		% [~, index_max] = max(likelihoods);
		% delta = h(index_max);

		Q.ProbWeights(end) = delta;
		Q.ProbWeights(index_min_support) = Q.ProbWeights(index_min_support) - delta;

		index_remove = (Q.ProbWeights == 0);

		Q.ProbWeights(index_remove) = [];
		Q.Support(index_remove) = [];
	end

	


end

function derivative = D(phi, X, Q, theta)
	f_Q = Q.ProbWeights * phi(X, Q.Support)';
	f_theta = phi(X,theta)';
	denom = repmat(f_Q, [length(theta), 1]);
	derivative = sum(f_theta ./ denom - 1, 2);
end

function [theta_star,derivative] = max_theta2(phi, X, Q)
	DQ = @(theta) D(phi, X, Q, theta);

	% Coarse
	res = 1000;
	xx = linspace(min(X), max(X), res);
	dx = xx(2) - xx(1);
	yy = DQ(xx);
	figure(2)
	plot(xx, yy)
	drawnow
	[~, ii] = max(yy);

	%Fine
	res2 = 1000;
	xx = linspace(xx(ii) - dx, xx(ii) + dx, res2);
	yy = DQ(xx);
	[derivative, ii] = max(yy);
	theta_star = xx(ii);
end

function log_likelihood = L(phi, X, Q)
	log_likelihood = sum(log(Q.ProbWeights * phi(X, Q.Support)'));
end

function likelihood = blah(phi, X, Q, index_min_support, x)
	Q_test = Q;
	Q_test.ProbWeights(end) = x;
	Q_test.ProbWeights(index_min_support) = Q_test.ProbWeights(index_min_support) - x;
	likelihood = L(phi, X, Q_test);
end