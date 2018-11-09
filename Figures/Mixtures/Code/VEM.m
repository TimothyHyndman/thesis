% X is column vector

function [Q,log_likelihood] = VEM(phi,X)

	%Initial value
	Q.Support = mean(X);
	Q.ProbWeights = 1;

	% m = 10;
	% Q.Support = linspace(min(X), max(X), m);
	% Q.ProbWeights = (1 / m) * ones(1, m);  

	looping = true;

	while looping
		% Find theta* = argmax D(theta)
		[theta_star, derivative] = max_theta2(phi, X, Q);
		theta_star

		if derivative < 1e-2
			looping = false;
		end

		% Add theta* to support
		Q.Support = [Q.Support, theta_star];
		Q.ProbWeights = [Q.ProbWeights, 0];

		% Find point of support theta- which minimizes D(theta_j)
		[~, index_min_support] = min(D(phi, X, Q, Q.Support));

		%Find how much mass to move from theta- to theta*
		res = 100;
		h_max = Q.ProbWeights(index_min_support);
		h = linspace(0, h_max, res);
		delta_h = h(2) - h(1);
		likelihoods = zeros(1, res);
		for i = 1:res
			likelihoods(i) = blah(phi, X, Q, index_min_support, h(i));
		end
		[~, index_max] = max(likelihoods);
		coarse_h_sol = h(index_max)

		% figure(1)
		% plot(h, likelihoods)
		% drawnow
		% pause

		res = 100;
		
		h = linspace(max(coarse_h_sol - delta_h, 0), min(coarse_h_sol + delta_h, h_max), res);
		likelihoods = zeros(1, res);
		for i = 1:res
			likelihoods(i) = blah(phi, X, Q, index_min_support, h(i));
		end
		[lll, index_max] = max(likelihoods);
		h_sol = h(index_max)

		% figure(1)
		% plot(h, likelihoods)
		% drawnow
		% pause

		

		% h_sol = Q.ProbWeights(index_min_support);


		Q.ProbWeights(end) = h_sol;
		Q.ProbWeights(index_min_support) = Q.ProbWeights(index_min_support) - h_sol;

		index_remove = (Q.ProbWeights == 0);

		Q.ProbWeights(index_remove) = [];
		Q.Support(index_remove) = [];
		Q.ProbWeights = Q.ProbWeights / sum(Q.ProbWeights);

		Q
		% lll
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
	hold on
	scatter(Q.Support, Q.ProbWeights)
	hold off
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