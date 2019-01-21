%---------------------------------------------------------------
%---------------------------------------------------------------
% Code written by Aurore Delaigle for the paper: Delaigle, A. and Hall, P. (to appear). Methodology for nonparametric deconvolution when the error distribution is unknown.  JRSSB  
% This is NOT the code used in the paper
% This is an attempt at a cleaned up version of the codes used in the paper, which might contain errors
% Do not distribute unless authorised by the author
% Contact Aurore Delaigle by email if you find errors in the code
%---------------------------------------------------------------
%---------------------------------------------------------------



function [c,ceq]=mycon(pval)

global fmax Term1boot Termequality Termequalitymax penalite_g penalite_g_max;
	
	%Termboot1 is a global variable computed within fobjBoot
	%Termequality is a global variable computed within fobjBoot

	%We want to solve the minimisation problem under the constraint that
	%Term1boot<=fmax:  the variance of X we are trying the minimise should be such that the integral (objective function) of the main problem is no larger than value found under var X constraint
	%Termequality<=0 (this=0 since it is an absolute value)
	
	ceq=[];
	c=[Term1boot-fmax, Termequality - Termequalitymax, penalite_g - penalite_g_max];

end
