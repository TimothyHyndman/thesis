%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------


function y=fobjBoot(pvec)
global fmax Term1boot fctargs Termequality penalite_g;

	fctarg=fctargs;
	xgrid=fctarg.xgrid;
	n=fctarg.n;
	tt=fctarg.tt;
	rehatphiW=fctarg.rehatphiW;
	imhatphiW=fctarg.imhatphiW;
	normhatphiW=fctarg.normhatphiW;
	sqabshatpsi=fctarg.sqabshatpsi;
	varW=fctarg.varW;

	pvec=reshape(pvec,1,length(pvec));
	pvec=[pvec,1-sum(pvec)];
	
	%Compute characteristic fucntion of discrete approximation
	
	OO=outerop(tt,xgrid,'*');
	pmat=repmat(pvec,length(tt),1);
	cosO=cos(OO).*pmat;
	sinO=sin(OO).*pmat;
	clear OO;
	
	rephip=sum(cosO,2);
	imphip=sum(sinO,2);

	%Since f_U is symemtric, then our estimator of phi_U needs to be real ==> Its imaginary part needs to be equal to zero.
	%Since we estimate phi_U by \hat phi_W/\phi_p, we need to have \hat phi_W/\phi_p is real, which we can replace by 
	%\hat phi_W \bar \phi_p/|\phi_p|^2 is real
	%or again \hat phi_W \bar \phi_p is real
	%or again imaginary part of \hat phi_W \bar \phi_p is equal to zero for all t. This is what I compute below:
	Termequality=sum(abs(rephip.*imhatphiW-imphip.*rehatphiW));

	normphip=sqrt(rephip.^2+imphip.^2);
	Realterm=rehatphiW.*normphip-sqabshatpsi.*rephip;
	Imterm=imhatphiW.*normphip-sqabshatpsi.*imphip;
    
    %impose a penalty if |phi_U| is greater than 1:
	hatphiU=normhatphiW./(normphip);
    penalite_g = sum(hatphiU(hatphiU>1) - 1);
	
	%Compute the variance of the discrete distribution
	EX=sum(xgrid.*pvec);
	hatvarX= sum((xgrid-EX).^2.*pvec);
	penalite2=hatvarX;
	
	%Compute value of the integral (objective function of the main optimisation problem)
	sigt=-tt(1)/2;
	weight=0.75/(2*sigt)*(1-(tt/(2*sigt)).^2);
	dt=tt(2)-tt(1);
	Term1boot=sum((Realterm.^2+Imterm.^2).*weight')*dt;
	y=penalite2;

