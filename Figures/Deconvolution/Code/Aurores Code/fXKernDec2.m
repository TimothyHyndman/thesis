%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------



function fXdec=fXKernDec2(xx,hPI,W,tlim,ppphiU,hatvarU)

	phiK = @(t) (1-t.^2).^3;
	muK2 = 6;

	dt = .0002;
	t = (-1:dt:1)';
	t=reshape(t,length(t),1);
	
	phiU=phiUspline(t/hPI,hatvarU,tlim,ppphiU);
	OO=outerop(t/hPI,W,'*');
	%Estimate empirical characersitic fucntion of W
	
	n=length(W);
	rehatphiW=sum(cos(OO),2)/n;
	imhatphiW=sum(sin(OO),2)/n;

	rephip=rehatphiW./phiU;
	imphip=imhatphiW./phiU;

	xt=outerop(t/hPI,xx,'*');
	fXdec=cos(xt).*repmat(rephip,1,length(xx))+sin(xt).*repmat(imphip,1,length(xx));
	fXdec=sum(fXdec.*repmat(phiK(t),1,length(xx)),1)/(2*pi)*dt/hPI;
