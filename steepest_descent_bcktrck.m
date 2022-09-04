function [xk, fk, gradfk_norm, k, xseq, btseq] = steepest_descent_bcktrck ...
(x0, f, gradf, alpha0, kmax, tolgrad, c1, rho, btmax, mem)
% Function that performs the steepest descent optimization method with
% backtracking strategy for the line search
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% alpha0 = starting value of alpha for backtracking
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
% mem = if equal to 1 then the method stores all xk values in xseq
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = array with the value f(xk) for each iteration k;
% gradfk_norm = values of the norm of gradf for each k
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
    if nargin == 9
        mem = 0;
    end    
    % Function handle for the armijo condition
    farmijo = @(fk, alpha, gradfk, pk) fk + c1 * alpha * gradfk' * pk;    
    % Initializations
    if mem == 1
        xseq = zeros(length(x0), kmax);
    else
        xseq = 0;
    end
    btseq = zeros(kmax, 1);
    xk = x0;
    fk = zeros(1, kmax+1);
    fk(1) = f(xk);
    gradfk = gradf(xk);
    k = 0;
    gradfk_norm = zeros(1, kmax+1);
    gradfk_norm(1) = norm(gradfk);
    while k < kmax && gradfk_norm(k+1) >= tolgrad
        % Compute the descent direction
        pk = -gradfk;   
        alpha = alpha0;
        % Compute the new value for xk
        xknew = xk + alpha * pk;
        fknew = f(xknew);    
        % backtracking for line search
        bt = 0;
        while bt < btmax && fknew > farmijo(fk(k+1), alpha, gradfk, pk)
            alpha = rho * alpha;
            xknew = xk + alpha * pk;
            fknew = f(xknew);
            bt = bt + 1;
        end
        % updating values
        xk = xknew;
        gradfk = gradf(xk);
        % Increase the step by one
        k = k + 1;
        % Compute the gradient of f in xk
        gradfk_norm(k+1) = norm(gradfk);
        fk(k+1) = fknew;  
        if mem == 1
            % Store current xk in xseq
            xseq(:, k) = xk;
        end
        btseq(k) = bt;
    end
    % Cut arrays to the correct size
    if mem == 1
        xseq = xseq(:, 1:k);
    end
    btseq = btseq(1:k);
    gradfk_norm = gradfk_norm(1:k+1);
    fk = fk(1:k+1);
end