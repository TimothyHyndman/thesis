%% n = size of sample, NSR = noise signal ratio

function [W,truedens,X,U, truepmf] = generatedata(n,NSR,dist_type,error_type,seed)

%Initiate random seeds to be able to replicate the results
if nargin == 5
    rng(seed);
end

truepmf = [];

%Generate sample of size n from X
switch dist_type
    case 'chi'
        dchi=3;
        X=chi2rnd(dchi,1,n);
        varX=2*dchi;
        X=X/sqrt(varX);
        varX=1;
        varXlong=2*dchi;
        cvar=sqrt(varXlong);
        truedens=@(xx) cvar*chi2pdf(xx*cvar,dchi);
    case 'gamma'
        shape_par = 5;
        scale_par = 2;
        X = gamrnd(shape_par,scale_par,[1,n]);
        varX = shape_par*scale_par^2;
        X=X/sqrt(varX);
        varX=1;
        varXlong = shape_par*scale_par^2;
        cvar=sqrt(varXlong);
        truedens=@(xx) cvar*gampdf(xx*cvar,shape_par,scale_par);
    case 'pareto'
        k=-0.5;
        varX=1;
        theta=2;
        X = gprnd(k,varX,theta,[1,n]);
        varX = varX^2/((1-k)^2*(1-2*k));
        X=X/sqrt(varX);
        varX=1;
        varXlong = varX^2/((1-k)^2*(1-2*k));
        cvar=sqrt(varXlong);
        truedens=@(xx) cvar*gppdf(xx*cvar,k,varX,theta);
    case 'normal'   %This is symmmetric so should stuff things up
        mu = 2;
        varX = 1;
        X=normrnd(mu,varX,[1,n]);
        varX=varX^2;
        X=X/sqrt(varX);
        varX=1;
        varXlong = varX^2;
        cvar=sqrt(varXlong);
        truedens=@(xx) cvar*normpdf(xx*cvar,mu,varX);
    case 'degenerate'
        a = 2;
        X = zeros(1,n);
        X = X+a;
        varX = 1;   %This isn't true!! But we can't use varX = 0 for the error bit
        truedens = @(xx) abs(xx - a) < 0.05;
    case 'discrete'
        X = zeros(1,n);
        for i = 1:n
            fred = unifrnd(0,1);
            if fred < 0.3
                X(i) = 1;
            else if fred < 0.5
                    X(i) = 2;
                else if fred < 0.6
                        X(i) = 3;
                    else if fred < 0.7
                            X(i) = 4;
                        else if fred < 0.75
                                X(i) = 6;
                            else if fred < 0.8
                                    X(i) = 7;
                                else if fred < 0.95
                                        X(i) = 9;
                                       else X(i) = 9.5;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        scale = sqrt(var(X));
%         X = X/scale;
        varX = 1;
        truedens = @(xx) abs(xx - 1/scale) < 0.05 | abs(xx - 2/scale) < 0.05 | abs(xx - 3/scale) < 0.05;
    case 'discrete2'
        X = zeros(1,n);
        for i = 1:n
            fred = unifrnd(0,1);
            if fred < 0.3
                X(i) = 1;
            else if fred < 0.5244
                    X(i) = 2.01;
                else if fred < 0.6
                        X(i) = 3.1416;
                    else
                        X(i) = 4;
                    end
                end
            end
        end
        scale = sqrt(var(X));
        scale = 1;
        X = X/scale;
        varX = 1;
        truedens = @(xx) abs(xx - 1/scale) < 0.05 | abs(xx - 2.01/scale) < 0.05 | abs(xx - 3.1416/scale) < 0.05 | abs(xx - 4/scale) < 0.05;
        truepmf.Support = [1, 2.01, 3.1416, 4] / scale;
        truepmf.ProbWeights = [0.3, 0.5244-0.3, 0.6-0.5244, 1- 0.6];
    case 'gumbel'
        mu = -1;
        varX = 1;
        sigma = sqrt(6*varX/pi^2);
        X = -evrnd(mu, sigma, [1,n]);
        truedens = @(xx) evpdf(-xx, mu, sigma);
end

%Choose error type (lap or norm) and compute true phi_U and var(U) for this error
switch error_type
    case 'norm'
        sigU=sqrt(NSR*varX);
        U=normrnd(0,sigU,1,n);
    case 'lap'
        sigU=sqrt(NSR*varX/2);
        U=rlap(sigU,1,n);
    case 'cauchy'
        a = sqrt(NSR*varX)/10;
        U = a.*randn(1,n)./randn(1,n);
    case 'discrete'
        U = zeros(1,n);
        varU = NSR*varX;
        for i = 1:n
            fred = unifrnd(0,1);
            if fred < 0.6
                U(i) = 0;
            elseif fred < 0.8
                U(i) = -varU/0.4;
            else
                U(i) = varU/0.4;
            end
        end
    case 'degenerate'
        U = zeros(1,n);
end

%Generate contaminated data.
W=X+U;