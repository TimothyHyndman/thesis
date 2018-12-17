%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------



function hPI = PI_deconvUestth4(W,tlim,ppphiU,hatvarU)

% Input variables:
%       W  - argument (n x 1 vector)
%       para  - a Matlab structure of input parameters
%               kernel = type of kernel,
%                       1: Gaussian kernel 
%                       2: phiK(t) = (1-t^2)^3 for -1<t<1  (default)
%               errdist = distribution of the error variable
%                       1: Normal (0,sigma^2)
%                       2: Laplace(0, lambda) 
%               sigma   = standard deviation of the Laplacel error
%               hgrid   = candidates for the optimal bandwidth 
%               th4     = theta4 in Delaigle & Gijbels (2004),
%                           \int (f^(4)(x))^2 dx, default = NR
% Output:
%       hPI = PI bandwidth which minimizes AMISE of the density estimator, see Delaigle and Gijbels
%
%Default values

phiK = @(t) (1-t.^2).^3;
muK2 = 6;
n = length(W);

%grid of h values where to search for a solution
maxh=(max(W)-min(W))/10;
hnaive=1.06*sqrt(var(W))*n^(-1/5);
hgrid=hnaive/3:(maxh-hnaive/3)/100:maxh;

lh = length(hgrid);

stdevx = max(sqrt(var(W) - hatvarU),1/n); % std(X)
th4 = stdevx^(-9)*105/(32*sqrt(pi)); % Estimate theta4 by NR     


dt = .0002;
t = (-1:dt:1)';
t=reshape(t,length(t),1);
hgrid=reshape(hgrid,1,lh);

toverh=t*(1./hgrid);

phiK2=(phiK(t)).^2;
phiU2=phiUspline(toverh,hatvarU,tlim,ppphiU).^2;


rr=3;
% Find h3 for th3
term1= -hgrid.^2*muK2*th4;
term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
term2=sum(term2,1)*dt./(2*pi*n*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;
indh3=find(ABias2==min(ABias2),1,'first');
h3 = hgrid(indh3);



OO=outerop(t/h3,W,'*');
%Estimate empirical characersitic fucntion of W
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
clear OO;
normhatphiW2=rehatphiW.^2+imhatphiW.^2;
th3 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh3))*dt/(2*pi*h3^(2*rr+1));




rr=2;
% Find h2 for th2
term1= -hgrid.^2*muK2*th3;
term2=repmat(t.^(2*rr).*phiK2,1,lh)./phiU2;
term2=sum(term2,1)*dt./(2*pi*n*hgrid.^(2*rr+1));

ABias2 = (term1 + term2).^2;
indh2=find(ABias2==min(ABias2),1,'first');
h2 = hgrid(indh2);


OO=outerop(t/h2,W,'*');
%Estimate empirical characersitic fucntion of W
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
clear OO;
normhatphiW2=rehatphiW.^2+imhatphiW.^2;
th2 = sum(t.^(2*rr) .* normhatphiW2 .* phiK2 ./ phiU2(:,indh2))*dt/(2*pi*h2^(2*rr+1));


term1=hgrid.^4*muK2^2*th2/4;
term2=repmat(phiK2,1,lh)./phiU2;
term2=sum(term2,1)*dt./(2*pi*n*hgrid);
AMISE=term1+term2;

indh=find(AMISE==min(AMISE),1,'first');
hPI = hgrid(indh);



end