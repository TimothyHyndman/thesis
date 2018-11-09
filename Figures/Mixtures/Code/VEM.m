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

		if derivative < 1e-2
			looping = false;
		end

		% Add theta* to support
		Q.Support = [Q.Support, theta_star];
		Q.ProbWeights = [Q.ProbWeights, 0];

		% Find point of support theta- which minimizes D(theta_j)
		[~, index_min_support] = min(D(phi, X, Q, Q.Support));
		index_reduce = index_min_support;
		index_increase = length(Q.Support);

		[Q, likelihood] = exchange_mass(phi, X, Q, index_reduce, index_increase);
		likelihood

		% Movement phase
		move_loop = true;
		while move_loop
			[Q.Support, I] = sort(Q.Support);
			Q.ProbWeights = Q.ProbWeights(I);

			support_diff = Q.Support(2:end) - Q.Support(1:end-1);
			[~, I] = min(support_diff);
			[~, which_index] = max(D(phi, X, Q, Q.Support(I:I+1)));

			if which_index == 2
				index_increase = I + 1;
				index_reduce = I;
			else
				index_increase = I;
				index_reduce = I + 1;
			end
			
			[Q_new, likelihood] = exchange_mass(phi, X, Q, index_reduce, index_increase);

			if length(Q_new.ProbWeights) == length(Q.ProbWeights)
				move_loop = false;
			else
				Q = Q_new;
			end
		end

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

function likelihood = blah(phi, X, Q, index_reduce, index_increase, x)
	Q_test = Q;
	Q_test.ProbWeights(index_increase) = Q_test.ProbWeights(index_increase) + x;
	Q_test.ProbWeights(index_reduce) = Q_test.ProbWeights(index_reduce) - x;
	likelihood = L(phi, X, Q_test);
end

function [Q, likelihood] = exchange_mass(phi, X, Q, index_reduce, index_increase)
	%Find how much mass to move from theta- to theta*
	res = 100;
	h_max = Q.ProbWeights(index_reduce);
	h = linspace(0, h_max, res);
	delta_h = h(2) - h(1);
	likelihoods = zeros(1, res);
	for i = 1:res
		likelihoods(i) = blah(phi, X, Q, index_reduce, index_increase, h(i));
	end
	[likelihood, index_max] = max(likelihoods);
	coarse_h_sol = h(index_max);

	% res = 100;
	% h = linspace(max(coarse_h_sol - delta_h, 0), min(coarse_h_sol + delta_h, h_max), res);
	% likelihoods = zeros(1, res);
	% for i = 1:res
	% 	likelihoods(i) = blah(phi, X, Q, index_reduce, index_increase, h(i));
	% end
	% [likelihood, index_max] = max(likelihoods);
	% h_sol = h(index_max);

	h_sol = coarse_h_sol;
	% figure(1)
	% plot(h, likelihoods)
	% drawnow
	% pause

	% Add that much mass to theta*
	Q.ProbWeights(index_increase) = Q.ProbWeights(index_increase) + h_sol;
	Q.ProbWeights(index_reduce) = Q.ProbWeights(index_reduce) - h_sol;

	index_remove = (Q.ProbWeights == 0);
	Q.ProbWeights(index_remove) = [];
	Q.Support(index_remove) = [];
	Q.ProbWeights = Q.ProbWeights / sum(Q.ProbWeights);
end