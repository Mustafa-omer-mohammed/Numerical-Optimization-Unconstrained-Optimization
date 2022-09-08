function [xk, fk, gradfk_norm, k, xseq, btseq] = ...
    newton_general(x0, f, gradf, Hessf, kmax, ...
    tolgrad, c1, rho, btmax, FDgrad, FDHess, h)
%
%
% [xk, fk, gradfk_norm, k, xseq] = ...
%     newton_general(x0, f, gradf, Hessf, ...
%     kmax, tolgrad, c1, rho, btmax, FDgrad, FDhess, h)
%
% Function that performs the newton optimization method, 
% implementing the backtracking strategy and, optionally, finite
% differences approximations for the gradient and/or the Hessian.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f (not necessarily
% used);
% Hessf = function handle that describes the Hessian of f (not necessarily 
% used);
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy;
% FDgrad = 'fw' (FD Forward approx. for gradf), 'c' (FD Centered approx. 
% for gradf), any other string (usage of input Hessf)
% FDHess = 'c' (FD Centered approx. for Hessf), any 
% other string (usage of input Hessf);
% h = approximation step for FD (if used).
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%

switch FDgrad
    case 'fw'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'fw')
        gradf = @(x) findiff_grad(f, x, h, 'fw');
        
    case 'c'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'c')        
        gradf = @(x) findiff_grad(f, x, h, 'c');
        
    otherwise
        % WE USE THE INPUT FUNCTION HANDLE gradf...
        %
        % THEN WE DO NOT NEED TO WRITE ANYTHING!
        % ACTUALLY WE COULD DELETE THE OTHERWISE BLOCK
        
end

switch FDHess
    case 'c'
        % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
        Hessf = @(x) findiff_Hess(f, x, sqrt(h));

end

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);

xk = x0;
fk = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);

while k < kmax && gradfk_norm >= tolgrad
    % Compute the descent direction as solution of
    % Hessf(xk) p = - graf(xk)
    pk = -Hessf(xk)\gradfk;
    
    
    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, xk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end
