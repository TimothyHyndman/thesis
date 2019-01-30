%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------


function reponse=rlap(szC,n1,n2)
% This function generates a number from a Laplace(szC)
y=unifrnd(0,1,n1,n2);
reponse= szC*log(2*y);
reponse(y>0.5)=-szC*log(2-2*y(y>0.5));
reponse;

