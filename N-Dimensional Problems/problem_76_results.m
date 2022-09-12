% Script for minimizing the function in problem 76
% INITIALIZATION
close all; clear; clc;
disp('** PROBLEM 76 **');
rho = 0.5; c = 1e-4; kmax = 10000; tolgrad = 1e-8;
alpha0 = 1;
btmax = 50;


% in_exact newton inputs
pcg_maxit = 50 ;
h = 1e-8 ;
FDgrad = '' ;  % we will use the exact gradient 
FDHess = '';   % we will use the exact hessian

% Function handles
f = @(x) problem_76_function(x);
gradf = @(x) problem_76_grad(x);
Hessf = @(x) problem_76_hess(x);

%------------------------------------------------------------------------%
% forcing terms 
select_fterm = ['l' , 's', 'q'] ;

% The Space Dimension 

n_values = [1e3,1e4,1e5]; 

%------------------------------------------------------------------------%

% 
for j = 1:length(n_values)
    n = n_values(j);
    methods = strings([5,1]);
    elapsed_times = zeros(5,1) ;
    grad_norm_last = zeros(5,1) ;
    fk_last = zeros(5,1) ;
    k_iterations = zeros(5,1) ;
    PCG_IT = zeros(3,1) ;
    format long
    
    disp(['SPACE DIMENSION: ' num2str(n, '%.0e')]);

    % generating starting point
    x0 = 2 * ones(n, 1);

    r = 0 ;
    value = j ;
    nor = j+3;
    steep = j+4;
    
    figure(value)  % create figure 1 Function value trend w.r.t num iteration k 
    figure(nor)  % create figure 1 norm gradient trend w.r.t num iteration k 
    colors = ['r',  'g' , 'y']; % colors for the plots

    for term = 1:length(select_fterm)
        ft_selected = select_fterm(term);
        fprintf('********************* fselected = %s ***************', ft_selected);
        switch ft_selected
            case 'l'  % linear convergence
                fterms =  @(gradfk,k)0.5 ;
                ft_print = ' IN-Newton linear' ;
            case  'q' % quadratic convergence
            fterms =  @(gradfk,k)min(0.5,norm(gradfk)) ;
            ft_print = ' IN-Newton quadratic' ;
            case  's' % super linear convergence
                fterms = @(gradfk,k)min(0.5,sqrt(norm(gradfk)));
                ft_print = ' IN-Newton super-inear' ;
    
            otherwise  % linear convergence by default
                fterms = @(gradfk,k)0.5 ;
                ft_print = ' Inexact Newton with linear rate of convergance' ;
        end
%         fprintf('********************* r = %f *************** \n', term);
        fprintf('*********************  %s ***************', ft_print);
            tic;
                [~, fk, gradfk_norm, k, ~, ~, pcg_iter,fk_seq] = ...
                innewton_general(x0, f, gradf, Hessf, kmax, ...
                tolgrad, c, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit) ;
                
            elapsed_time = toc;

        %--------------- Collecting results 

                methods(term,1) = ft_print ;% 1st column method name
                elapsed_times(term) = elapsed_time ; % second column elapsed time
                grad_norm_last(term) = gradfk_norm(end) ; % 3rd column grad_norm last
                fk_last(term) = fk(end) ; % 4th column function value last
                k_iterations(term) = k ; % 5th column k iteration
                PCG_IT(term) = sum(pcg_iter) ;
                r = term+1 ;
                figure(value),yyaxis left,semilogy(fk_seq, 'LineWidth', 2,'Color',colors(term),'DisplayName', ft_print ),hold on;
                figure(nor),yyaxis left,semilogy(gradfk_norm, 'LineWidth', 2,'Color',colors(term),'DisplayName', ft_print ),hold on;
    end
    M = max(k_iterations)+3;  % setting the X axis tickes 

    tic;
    [~, fk, gradfk_norm, k, ~, ~, gmres_it] = ...
        newton_bcktrck(x0, f, gradf, Hessf, kmax, tolgrad, c, rho, btmax);
    elapsed_time = toc;

    %--------------- Collecting results 

    methods(r) = 'Newton exact' ;% 1st column method name
    elapsed_times(r) = elapsed_time ; % second column elapsed time
    grad_norm_last(r) = gradfk_norm(end) ; % 3rd column grad_norm last
    fk_last(r) = fk(end) ; % 4th column function value last
    k_iterations(r) = k ; % 5th column k iteration
    r = r+1 ;

    %------------------- Plotting the grapgh for function value

    figure(value),yyaxis left,semilogy(fk, 'LineWidth', 2,'Color','b','DisplayName', 'Newton Exact'),hold off ,yyaxis left,
    ylabel('Value of the function'),
    xlabel('K iteration'),
    grid on 
    xticks(1:M),
    legend('Location', 'southwest') ,title(['Function value Trend comparing all methods for ', num2str(n) , ' Dimensions'] ,'FontSize',10 );
    temp1 = ['./76_images/Function-value',num2str(n),'.png'];

    %------------------- Exporting the grapgh for function value

    exportgraphics(figure(value),temp1,'Resolution',500) ;

    %------------------- Plotting the grapgh for Norm Gradient
    
    figure(nor),yyaxis left,semilogy(gradfk_norm, 'LineWidth', 2,'Color','b','DisplayName', 'Newton Exact'), grid on;
    hold off ,
    ylabel('Norm-Gradient'),
    xlabel('K iteration'),
    grid on 
    xticks(1:M),
    legend('Location', 'southwest') ,
    figure(nor),title(['Norm Gradient Trend comparing all methods for ', num2str(n) , ' Dimensions'] ,'FontSize',10);
    temp2 = ['./76_images/gard-norm',num2str(n),'.png'];

    %------------------- Exporting the grapgh for Norm Gradient

    exportgraphics(figure(nor),temp2,'Resolution',500) ;
    
    tic;
    [~, fk, gradfk_norm, k, ~, btseq] = ...
        steepest_descent_bcktrck(x0, f, gradf, alpha0, kmax, ...
            tolgrad, c, rho, btmax);
    elapsed_time = toc;

    %--------------- Collecting results 

    methods(r,1) = 'Steepest Decent' ;% 1st column method name
    elapsed_times(r) = elapsed_time ; % second column elapsed time
    grad_norm_last(r) = gradfk_norm(end) ; % 3rd column grad_norm last
    fk_last(r) = fk(end) ; % 4th column function value last
    k_iterations(r) = k ; % 5th column k iteration

    %----- Plotting the grapgh for Norm Gradient & Function value

    figure(steep),
    yyaxis left;
    semilogy(fk, 'LineWidth', 2);   % plot function value

    hold on 

    semilogy(gradfk_norm, '--', 'LineWidth', 1 , 'Color' , 'r'); % plot Norm grad
    ylabel('Norm of the gradient');
    legend('Fk trend', 'GradFk trend', 'Location', 'northeast');
    grid on 
    hold off
    title(['Function value Trend Steepest Decent ', num2str(n) , ' Dimensions'] ,'FontSize',10 );
    temp3 = ['./76_images/Sttepest-Decent',num2str(n),'.png'];
    %----- Exproting the grapgh for Norm Gradient & Function value

    exportgraphics(figure(steep),temp3,'Resolution',500) ;

    % --------------------------------- printing Resutls ---------------------%

    disp(['Results with N =  ', num2str(n) , ' Dimensions ']);
    disp(['% ---------------------------------','Results with N =  ', num2str(n) , ' Dimensions' ,' ---------------------%'])
    disp(['            Method             ','k          ','elapsed_times   ', 'grad_norm_last    ' , 'fk_last      ' ]),
    disp([methods ,  fix(k_iterations), double(elapsed_times),double(grad_norm_last) , double(fk_last)]) ;
    
end  

