function [Q,log_likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,X)
% Mixing Homogeneity Paper Log Likelihood Stuff

%Example usage: 
% sigma = 0.5;
% phi = @(x,u) (1/(sigma*sqrt(2*pi)))*exp(-((repmat(x,[1,length(u)]) - repmat(u,[length(x),1])).^2)./(2*sigma^2));
% Y = unifrnd(0,1,[50,1]);
% [Q, log_likelihood, Q_unsimplified] = MixtureLikelihoodMovingMasses2(phi,Y);

m = min(10,length(X));
a = min(X);
b = max(X);

Q.Support = linspace(a,b,m);

if m == 1
    Q.Support = mean(X);
elseif length(X) <= 3
    Q.Support = X';
    m = length(X);
end

Q.ProbWeights = (1/m)*ones(1,m);  

x0 = [Q.Support,Q.ProbWeights]; 
func = @(x) -L(phi,X,[x,1 - sum(x((end+1)/2 + 1:end))]);
x0(end) = [];
log_likelihood = -func(x0);

%Make sure that initial point is defined
if func(x0) == Inf
    display('Initial point undefined')
    looping = true;
    while looping
        Q.Support = unifrnd(a,b,[1,m]);
        x0 = [Q.Support,Q.ProbWeights];
        x0(end) = [];
        display(x0)
        if func(x0) < Inf
            looping = false;
            log_likelihood = -func(x0);
        end
    end
end

options = optimoptions('fmincon','Algorithm','sqp','Display','off','TolFun',1.0e-09,'TolX',1.0e-09);  %sqp algorithm doesn't get stuck on some problems that default does
tolerance = 1.0e-06;   %Any lower and it gets stuck
counter = 1;
looping = true;

while looping
    %Summation and boundary constraints
    m = (length(x0)+1)/2;
    lb_points = a*ones([1,m]) - 0.01;
    ub_points = b*ones([1,m]) + 0.01;
    A_points = zeros([1,m]);
    
    lb_masses = zeros([1,m-1]);
    ub_masses = ones([1,m-1]);
    A_masses = ones([1,m-1]);

    lb = [lb_points,lb_masses];
    ub = [ub_points,ub_masses];
    A = [A_points,A_masses];
    b = 1;

    %Solve
    [solution_x,fval] = fmincon(func,x0,A,b,[],[],lb,ub,[],options);
    if -fval >= log_likelihood
        log_likelihood = -fval;
        Q.Support = solution_x(1:m);
        Q.ProbWeights = [solution_x(m+1:end),1 - sum(solution_x(m+1:end))];
        Q_unsimplified = Q;
        Q = SimplifyDist(Q);
    else
%         display(['fmincon ',num2str(counter),' failed'])
    end

    [theta_star,derivative] = max_theta2(phi,X,Q);

    if derivative <= 0 + tolerance
        break
    else
        [Q_test,log_likelihood_test] = max_weights(phi,X,Q,theta_star);
        Q_unsimplified = Q;
        Q_test = SimplifyDist(Q_test);
        if log_likelihood_test >= log_likelihood
            Q = Q_test;
            log_likelihood = log_likelihood_test;
            x0 = [Q.Support,Q.ProbWeights];
            x0(end) = [];
            Q_unsimplified = Q;
            Q = SimplifyDist(Q);
            [~,derivative] = max_theta2(phi,X,Q);
            if derivative <= 0+ tolerance
                break
            end
        else
%             display(['max weights ',num2str(counter),' failed'])
        end
    end

%Used for testing points that won't converge
    counter = counter+1;
    if counter > 10
%         display(X)
        display('FAILED')  
%         [theta_star,derivative] = max_theta2(phi,X,Q);
%         display([theta_star,derivative])
        break
    end
end

% plotDD(phi, X, Q)
end

%Subfunctions
function log_likelihood = L(phi,X,x)
Support = x(1:length(x)/2);
ProbWeights = x(length(x)/2+1:end);
log_likelihood = sum(log(ProbWeights*phi(X,Support)'));
end

function derivative = D(phi,X,theta,Q)
f_Q = Q.ProbWeights*phi(X,Q.Support)';
f_theta = phi(X,theta)';
denom = repmat(f_Q,[length(theta),1]);
derivative = sum(f_theta./denom - 1,2);
end

function [theta_star,derivative] = max_theta(phi,X,Q)
DQ = @(theta) D(phi,X,theta,Q);
res = 100000; %Need this high res near solution
xx = linspace(min(X),max(X),res);
yy = DQ(xx);
[derivative,ii] = max(yy);
theta_star = xx(ii);
end

function [theta_star,derivative] = max_theta2(phi,X,Q)
DQ = @(theta) D(phi,X,theta,Q);
res = 1000;
% As 
xx = linspace(min(X),max(X),res);
dx = xx(2) - xx(1);
yy = DQ(xx);
[~,ii] = max(yy);
%Fine
res2 = 1000;
xx = linspace(xx(ii) - dx,xx(ii)+dx,res2);
yy = DQ(xx);
[derivative,ii] = max(yy);
theta_star = xx(ii);
end

function [Q,log_likelihood] = max_weights(phi,X,Q,theta_star)
func = @(weights) -L(phi,X,[Q.Support,theta_star,weights,1 - sum(weights)]);   %function to minimize
lb = zeros([1,length(Q.Support)]);
A = ones([1,length(Q.Support)]);
b = 1;
options = optimoptions('fmincon','Algorithm','sqp','Display','off','TolFun',1.0e-09,'TolX',1.0e-09);
x0 = Q.ProbWeights;
% lh1 = -func(x0)
[probweights,fval] = fmincon(func,x0,A,b,[],[],lb,[],[],options);
% 1 - sum(probweights)
% lh2 = -func(probweights)
log_likelihood = -fval;
Q.Support = [Q.Support,theta_star];
Q.ProbWeights = [probweights,1 - sum(probweights)];
end

function plotDD(phi,X,Q)
DQ = @(theta) D(phi,X,theta,Q);
res = 1000;
xx = linspace(min(X),max(X),res);
yy = DQ(xx);
plot(xx,yy)
end

function Q = SimplifyDist(Q)
grid_points = Q.Support;
grid_masses = Q.ProbWeights;
zero_tolerance = 0.0000005;  %Need this at least as small as 0.0000005
search_size = 0.05;  %Getting points close to each other is rare - unlike getting points with low mass.
[grid_points,I] = sort(grid_points);
grid_masses = grid_masses(I);

%Get rid of zeros
index = 1;
looping = true;
while looping
    if grid_masses(index) < zero_tolerance
        grid_points(index:end-1) = grid_points(index+1:end);
        grid_points(end) = [];
        grid_masses(index:end-1) = grid_masses(index+1:end);
        grid_masses(end) = [];
    else
        index = index+1;
    end
    if index >= length(grid_points)
        looping = false;
    end
end

if grid_masses(end) < zero_tolerance
    grid_masses(end) = [];
    grid_points(end) = [];
end

%Get rid of duplicates
index = 1;
looping = true;
if length(grid_points) < 2
    looping = false;
end
while looping
    if abs(grid_points(index) - grid_points(index+1)) < search_size
        grid_points(index) = (grid_points(index)*grid_masses(index)+grid_points(index+1)*grid_masses(index+1))/(grid_masses(index)+grid_masses(index+1));  %Weighted average location
        grid_points(index+1:end-1) = grid_points(index+2:end);
        grid_masses(index) = grid_masses(index)+grid_masses(index+1);
        grid_masses(index+1:end-1) = grid_masses(index+2:end);
        grid_masses(end) = [];
        grid_points(end) = [];
    else
        index = index+1;
    end
    if index == length(grid_points)
        looping = false;
    end
end

grid_masses = grid_masses/sum(grid_masses);    %Normalise
Q.Support = grid_points;
Q.ProbWeights = grid_masses;
end