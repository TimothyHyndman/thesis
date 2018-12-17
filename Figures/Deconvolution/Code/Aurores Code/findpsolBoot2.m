%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------


function [psol,fval]=findpsolBoot2(m,objfunc)

global fctargs

options = optimoptions('fmincon','Display','off','Algorithm','active-set','TolFun',1e-6);  


%Setup constraints in form A*psol<b. We need that entries of x add up to 1 and
%are all non-negative.
A=-eye(m-1);
A=[A;ones(1,m-1)];
b=zeros(1,m-1);
b=[b,1];
b=b';

%Generic start where each point is given the same mass.
pstart=repmat(1/m,m,1);
pstart(m)=[];


[psol,fval]=fmincon(str2func(objfunc),pstart,A,b,[],[],[],[],[],options);

psol=[psol',1-sum(psol)];


%Find solution with 20 random starting points
for i=1:20

	pstart=psol.*unifrnd(1/2,2,[1,m]);  %Generate random start
	pstart=pstart./sum(pstart); %Normalize it
	pstart(m)=[];
	
	[psolnew,fvalnew]=fmincon(str2func(objfunc),pstart,A,b,[],[],[],[],[],options);

	if(fvalnew<fval)
	
		fval=fvalnew;
		psol=[psolnew,1-sum(psolnew)];
	end
end
