function [ x ] = asls( A, b, p, x0, nn, bounds)
% function [ x ] = asls( A, b, p )
% Asymmetric least squares.  Algorithm finds coefficients x of a linear
% model, A*x = b, where b is the measurement vector and A is a matrix of
% pure components.  Starting from an initial, least-squares estimate of x,
% the program finds a residual vector r = b-A*x. A merit function f is
% constructed, with f = sum(w.*r.^2), where w is a weighting vector that
% has values of p for positive residuals and 1-p for negative residuals.
% This asymmetric penalty function helps to ensure positive-going
% residuals.  p is a critical parameter and must be set by the user.  See,
% for example, Journal of Chromatography A, 1057 (2004) 21-30.
%
% INPUTS:   A   =   matrix of pure components.
%                   [m x n]
%           b   =   measurement vector.
%                   [m x 1]
%           p   =   penalty term (0 < p < 1)
%                   [1 x 1]
%           x0  =   initial guess of x0 (like from a standard LS fit
%                   [n x 1]
%           nn  =   use nonnegative ls for initial guess
%                   [1 x 1] boolean
%           bounds= limits to estimates of x (to prevent ridiculous values)
%                   [n x 2] (in the style of [lower, upper]
%
% OUTPUTS:  x   =   coefficient vector
%                   [n x 1]
%
% CHANGELOG
%   10/04/2012 - created function - zjs

if nargin<=4
    nn = 0 ;
end

% make an initial guess of x

if isempty(x0)
    if nn
        x0 = lsqnonneg(A,b);
    else
        x0 = A\b;
    end
end

% optimize, starting from the initial guess
if nargin<6
    bounds = -1;
end

options = optimset('MaxFunEvals',10000,'MaxIter',10000);
x = fminsearch(@(x) aslsres(x,A,b,p,bounds),x0,options);
    
    function f = aslsres(x,A,b,p,bounds)
        % f = merit function to evaluate fitness of candidate solution x.
        r = b-A*x;
        % compute weightings
        w(r>=0) = p;
        w(r<0) = 1-p;

        % calculate f
        f = sum(w'.*r.^2);
        if size(bounds,2)==2
            if any(x<bounds(:,1) | x>bounds(:,2));
                f=9e99;
                fprintf('9e99\n')
            end
        end

    end

end

