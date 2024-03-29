%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------

function [rephip,imphip,normphip]=computephiX(tt,xgrid,psol,n);


% % AURORE'S METHOD
% OO=outerop(tt,xgrid,'*');
% pmat=repmat(psol,length(tt),1);
% cosO=cos(OO).*pmat;
% sinO=sin(OO).*pmat;
% 
% old_rephip=sum(cosO,2);
% old_imphip=sum(sinO,2);

% TIM'S METHOD
OO = reshape(tt, [], 1) * reshape(xgrid, 1, []);
rephip = cos(OO) * psol(:);
imphip = sin(OO) * psol(:);



normphip=sqrt(rephip.^2+imphip.^2);
