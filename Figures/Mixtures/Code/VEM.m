% X is column vector

function [Q,log_likelihood] = VEM(phi, X, d_phi)

	%Initial value
	Q.Support = mean(X);
	Q.ProbWeights = 1;
	likelihood = Inf;
	% m = 10;
	% Q.Support = linspace(min(X), max(X), m);
	% Q.ProbWeights = (1 / m) * ones(1, m);  

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
		index_decrease = index_min_support;
		index_increase = length(Q.Support);

		% exchange_mass_2(phi, X, Q, index_decrease, index_increase);
		[Q, likelihood_new] = exchange_mass(phi, X, Q, index_decrease, index_increase);
		
		% pause
		% likelihood_new - likelihood
		% likelihood = likelihood_new;

		% % Movement phase
		% move_loop = true;
		% while move_loop
		% 	[Q.Support, I] = sort(Q.Support);
		% 	Q.ProbWeights = Q.ProbWeights(I);

		% 	support_diff = Q.Support(2:end) - Q.Support(1:end-1);
		% 	[~, I] = min(support_diff);
		% 	[~, which_index] = max(D(phi, X, Q, Q.Support(I:I+1)));

		% 	if which_index == 2
		% 		index_increase = I + 1;
		% 		index_decrease = I;
		% 	else
		% 		index_increase = I;
		% 		index_decrease = I + 1;
		% 	end
			
		% 	[Q_new, likelihood] = exchange_mass(phi, X, Q, index_decrease, index_increase);

		% 	if length(Q_new.ProbWeights) == length(Q.ProbWeights)
		% 		move_loop = false;
		% 	else
		% 		Q = Q_new;
		% 	end
		% end

		% Collection phase
		%IDEA: TRY TO COLLECT POINTS IN SAME MODE INTO ONE POINT
		move_loop = true;
		while move_loop
			[Q.Support, I] = sort(Q.Support);
			Q.ProbWeights = Q.ProbWeights(I);
			support_diff = Q.Support(2:end) - Q.Support(1:end-1);
			[~, I] = min(support_diff);


			iii = [I, I+1];
			Q_new = collate_points(phi, X, Q, iii);

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
	[derivative, ii] = max(yy);

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

function likelihood = likelihood_Q1_to_Q2(phi, X, Q1, Q2, pp)
	likelihood = sum(log(((1 - pp)' * Q1.ProbWeights + pp' * Q2.ProbWeights) * phi(X, Q1.Support)'), 2);
end


function [Q, likelihood] = exchange_mass(phi, X, Q, index_decrease, index_increase)
	%Find how much mass to move from theta- to theta*
	res = 1000;
	pp = linspace(0, 1, res);
	Q1 = Q;
	Q2 = Q;
	Q2.ProbWeights(index_increase) = Q2.ProbWeights(index_increase) + Q2.ProbWeights(index_decrease);
	Q2.ProbWeights(index_decrease) = 0;

	likelihoods = likelihood_Q1_to_Q2(phi, X, Q1, Q2, pp);

	[likelihood, index_max] = max(likelihoods);
	p_sol = pp(index_max);

	% figure(1)
	% plot(pp, likelihoods)
	% drawnow
	
	Q.ProbWeights = (1 - p_sol) * Q1.ProbWeights + p_sol * Q2.ProbWeights;

	index_remove = (Q.ProbWeights == 0);
	Q.ProbWeights(index_remove) = [];
	Q.Support(index_remove) = [];
	Q.ProbWeights = Q.ProbWeights / sum(Q.ProbWeights);
end

function [Q, likelihood] = exchange_mass_2(phi, X, Q, index_decrease, index_increase)
	% Along Q_1 to Q_2, the likelihood is a concave function. We can use this to find the minimium faster
	Q1 = Q;
	Q2 = Q;
	Q2.ProbWeights(index_increase) = Q2.ProbWeights(index_increase) + Q2.ProbWeights(index_decrease);
	Q2.ProbWeights(index_decrease) = 0;
	
	res = 100;
	pp = linspace(0, 1, res);

	f_Qpp = ((1 - pp)' * Q1.ProbWeights + pp' * Q2.ProbWeights) * phi(X, Q.Support)';
	f_Q2 = Q2.ProbWeights * phi(X, Q2.Support)';
	derivative = sum((repmat(f_Q2, [res, 1]) - f_Qpp) ./ f_Qpp, 2);

	figure(3);
	plot(pp, derivative);
	drawnow

	Q = 1;
	likelihood = 1;

end

function Q_new = collate_points(phi, X, Q, iiii)
	Q_new = Q;
	current_likelihood = L(phi, X, Q);
	p_collate = sum(Q.ProbWeights(iiii));
	theta_min = min(Q.Support(iiii));
	theta_max = max(Q.Support(iiii));

	res = 1000;
	theta_collate = linspace(theta_min, theta_max, res);

	% likelihood = sum(log(((1 - pp)' * Q1.ProbWeights + pp' * Q2.ProbWeights) * phi(X, Q1.Support)'), 2);

	ps = Q.ProbWeights;
	thetas = Q.Support;

	ps(iiii) = [];
	thetas(iiii) = [];

	likelihood = zeros(1, res);
	for i = 1:res
		likelihood(i) = sum(log( [ps, p_collate] * phi(X, [thetas, theta_collate(i)])'));
	end


	[d_likelihood, theta_index] = max(likelihood - current_likelihood);
	if d_likelihood > 0
		Q_new.Support = [thetas, theta_collate(theta_index)];
		Q_new.ProbWeights = [ps, p_collate];
	end

	% figure(4)
	% plot(theta_collate, likelihood - current_likelihood)
	% drawnow

	
end