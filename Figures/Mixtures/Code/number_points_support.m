function m = number_points_support(phi, X, Q_hat)
	%IDEA: Count the number of local maximums in D_Q(\theta) to determine number 
	% of points of support in Q_hat. I suspect this is more reliable than just
	% counting the non-zero weights.

	res = 1000;
	w = max(X) - min(X);

	theta = linspace(min(X) - 0.5*w - 0.1, max(X)+0.5*w+0.1, res);

	derivative = D(phi, X, Q_hat, theta);

	d_derivative = derivative(2:end) - derivative(1:end-1);

	sign_d_derivative = d_derivative >= 0;

	m = (sum(xor(sign_d_derivative(2:end), sign_d_derivative(1:end-1))) + 1)/2;
end

function derivative = D(phi, X, Q, theta)
	f_Q = Q.ProbWeights * phi(X, Q.Support)';
	f_theta = phi(X,theta)';
	denom = repmat(f_Q, [length(theta), 1]);
	derivative = sum(f_theta ./ denom - 1, 2);
end