%-------------  Newton WITH BACK TRACKING  ------------- %

function [xk, fk, gradfk_norm, k, xseq, btseq, gmres_it] = ...
newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c1, rho, btmax, mem)
% Function that performs the newton optimization method, 
% implementing the backtracking strategy.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
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
    btseq = zeros(1, kmax);
    gmres_it = zeros(1, kmax);
    xk = x0;
    fk = zeros(1, kmax+1);
    fk(1) = f(xk);
    k = 0;
    gradfk = gradf(xk);
    gradfk_norm = zeros(kmax+1, 1);
    gradfk_norm(1) = norm(gradfk);
    while k < kmax && gradfk_norm(k+1) >= tolgrad
        % Compute the descent direction as solution of
        % Hessf(xk) pk = -gradfx
        HessFx = Hessf(xk);
        [pk, flag, ~, iter] = gmres(HessFx, -gradfk);

        if flag ~= 0
            disp(flag)
            warning("gmres did not converge");
        end
        % Reset the value of alpha
        alpha = 1;
        % Compute the candidate new xk
        xnew = xk + alpha * pk;
        % Compute the value of f in the candidate new xk
        fnew = f(xnew);
        bt = 0;
        
        % Backtracking strategy: 
        % 2nd condition is the Armijo condition not satisfied
        while bt < btmax && fnew > farmijo(fk(k+1), alpha, gradfk, pk)
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
        gradfk = gradf(xk);
        % Increase the step by one
        k = k + 1;
        gradfk_norm(k+1) = norm(gradfk);
        fk(k+1) = fnew;    
        if mem == 1
            % Store current xk in xseq
            xseq(:, k) = xk;
        end
        % Store bt iterations in btseq
        btseq(k) = bt;
        gmres_it(k) = iter(2);
    end
    % "Cut" arrays to the correct size
    if mem == 1
        xseq = xseq(:, 1:k);
    end
    btseq = btseq(1:k);
    gradfk_norm = gradfk_norm(1:k+1);
    gmres_it = gmres_it(1:k);
    fk = fk(1:k+1);
end