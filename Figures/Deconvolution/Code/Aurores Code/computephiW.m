%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------

function [tt1,tt2,rehatphiW,imhatphiW,normhatphiW]=computephiW(tt,longt,W,n)


OO=outerop(tt,W,'*');

%-----------------------------------------------
%Estimate empirical characersitic fucntion of W
%-----------------------------------------------

rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
normhatphiW=sqrt(rehatphiW.^2+imhatphiW.^2);


%------------------------------------------------------------
%t^*: Keep only t-values such that |\hat_phi_W|< n^{-0.25}
%------------------------------------------------------------

tmp=tt(normhatphiW<n^(-0.25));
%Refine the interval of t-values accordingly
tt1=max(tmp(tmp<0));
tt2=min(tmp(tmp>0));
if isempty(tmp(tmp<0))
	tt1=min(tt);
end

if isempty(tmp(tmp>0))
	tt2=max(tt);
end

tt=tt1:(tt2-tt1)/longt:tt2;


%-----------------------------------------------------------------------
%Estimate empirical characersitic fucntion of W on that refined interval
%-----------------------------------------------------------------------

OO=outerop(tt,W,'*');
rehatphiW=sum(cos(OO),2)/n;
imhatphiW=sum(sin(OO),2)/n;
normhatphiW=sqrt(rehatphiW.^2+imhatphiW.^2);