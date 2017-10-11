function [z] = whittaker_smoother(groupy,lambda)
% function [z] = whittaker_smoother(groupy,lambda)
% 
% This is a function to implement the Whittaker smoother presented by
% Eilers in his paper "A Perfect Smoother," Analytical Chemistry, v.75 
% pp. 3631-3636, 2003.
%
% The smoother operates on the concept of penalized least squares, where
% the penalty function attempts to balance two competing goals of roughness
% and data fidelity.  One can think of fitting a spectrum with a smoothed
% version of itself and taking the residual.  Squaring and summing the
% residuals measures the fidelity, S.  Squaring and summing the first order
% differential of the smoothed fit (sum((z(i)-z(i+1)).^2)) gives an
% estimate of the fit's roughness, R.  A penalty function could be Q = S +
% lambda*R, where lambda is a free parameter (often called the Lagrange
% multiplier) chosen by the user and determining the amount of data
% smoothing.  In the paper it is shown that one can construct an operator
% that operates on a vector y to smooth it to a new smoothed function z.
% That operator is inv(I+lambda*D'*D), where I is the identity matrix, and
% D is the first order derivative of I.
%
% The paper goes into some detail in a more-or-less objective way to
% determine lambda - making this method a potentially parameterless fit.
% However, that also makes the smoother computationally costly.  Instead,
% the function here leaves the choice of lambda to the user.
%
% The chief advantages of this smoother over the traditional Savitzky-Golay
% smoother is that it intrinsically is capable of handling endpoints and
% large-stretches of missing data.  It also has computational advantages 
% when the window size of the SG smoother becomes large or the spectrum to 
% be smoothed becomes long.
%
% INPUTS:   groupy  -   The spectrum to be smoothed.
%                       [1 x m], [m x 1] or [m x n]
%           lambda  -   The Lagrange multiplier
%                       [1 x 1]
%
% OUTPUTS:  z       -   The smoothed spectrum.
%                       [1 x m], [m x 1] or [m x n], will follow y's shape.

% CREATED -zjs 09/08/2010
% UPDATE LOG
%       02/24/2011, zjs - updated to process multiple spectra at once.

if size(groupy,1) == 1 && size(groupy,2) ~= 1
    groupy = groupy';
    flag = 1;
else
    flag = 0;
end

for ij = 1:size(groupy,2)
    y = groupy(:,ij);
    m = length(y);
    E = speye(m);
    D = diff(E);
    C = chol(E+lambda*D'*D);
    z(:,ij) = C\(C'\y);
end

if flag
    z = z';
end
