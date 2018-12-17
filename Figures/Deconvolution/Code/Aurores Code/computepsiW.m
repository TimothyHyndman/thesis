%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------

function [rehatpsi,imhatpsi,sqabshatpsi]=computepsiW(tt,W,n)


WW=outerop(W,W,'-');
WW=reshape(WW,1,n^2);
indice=1;
for i = 2:n
	indice=[indice,(i-1)*n+i];
end
WW(indice)=[];

rehatpsi=0*tt;
imhatpsi=0*tt;

for i=1:length(tt)

	rehatpsi(i)=sum(cos(tt(i)*WW));
	imhatpsi(i)=sum(sin(tt(i)*WW));
end

rehatpsi=rehatpsi'/n/(n-1);
imhatpsi=imhatpsi'/n/(n-1);

sqabshatpsi=sqrt(sqrt(rehatpsi.^2+imhatpsi.^2));
