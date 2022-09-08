% INITIALIZATION
close all;
clear;
clc;

rho = 0.5;
c = 1e-4;
kmax = 1000;
tolgrad = 1e-12;
btmax = 50;
alpha0 = 0.2;
FDgrad = '';  % we used the exact gradient 
FDHess = 'c'; % we used the approximated hessian with finite difference 
pcg_maxit = 100 ;
h = sqrt(eps);
%------------------------------------------------------------------------%
% select the forcing term here
ft_selected = 'l';
ft_print = 'f_terms with linear rate of convergance' ;
switch ft_selected
    case 'l'  % linear convergence
        fterms =  @(gradfk,k)0.5 ;
        ft_print = 'fterms with linear rate of convergance' ;
    case  'qu' % quadratic convergence
    fterms =  @(gradfk,k)min(0.5,norm(gradfk)) ;
    ft_print = 'fterms with quadratic rate of convergance' ;
    case  's' % super linear convergence
        fterms = @(gradfk,k)min(0.5,sqrt(norm(gradfk)));
        ft_print = 'fterms with super linear rate of convergance' ;

    otherwise  % linear convergence by default
        fterms = @(gradfk,k)0.5 ;
end
%------------------------------------------------------------------------%
% ROSENBROCK FUNCTION
ros_func = @(x)100*(x(2,:)-x(1,:).^2).^2+(1-x(1,:)).^2;
ros_grad = @(x)[400*x(1,:).^3-400*x(1,:).*x(2,:)+2*x(1,:)-2; ...
    200*(x(2,:)-x(1,:).^2)];
ros_hess = @(f, x, h)findiff_Hess(f, x, h);

%x0 = [1.2; 1.2];
x0 = [-1.2; 1];
%x0 = [-5; -5];
disp("x0: " + mat2str(x0));

disp("****** TESTING ROSENBROCK FUNCTION *******");
disp("**** Inexact NEWTON METHOD WITH BACKTRACKING *****");
tic;
[xk_1, fk_1, gradfk_norm, k_1, xseq_1, btseq_1] = ...
    innewton_general(x0, ros_func, ros_grad, ros_hess, kmax, ...
    tolgrad, c, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit);
ex = toc;

disp('************** FINISHED ***************');
disp('************** RESULTS ****************');
fprintf('Results for x0 [%.2f , %.2f] \n', x0(1),x0(2));
disp(['xk: ', mat2str(xk_1), ' (actual minimum: [1; 1]);']);
disp(['f(xk): ', num2str(fk_1), ' (actual min. value: 0);']);
disp(['N. of Iterations: ', num2str(k_1),'/',num2str(kmax), ';']);
fprintf('EXCUTION TIME IS : %.2f ms \n' , ex *100 );
disp('************************************');

t = 'Inexact NEWTON METHOD WITH BACKTRACKING %s' ;
%title(sprintf('y(x) vs x using %s', filename))
% PLOTS
% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));
Z = 100*(Y-X.^2).^2+(1-X).^2;

fig1_n = figure();
% Contour plot with curve levels for each point in xseq
[C, ~] = contour(X, Y, Z);
hold on;
% plot of the points in xseq
plot([x0(1) xseq_1(1, :)], [x0(2) xseq_1(2, :)], '--*');
hold off;
title(sprintf(t, ft_print));
%title('Newton with backtracking - Rosenbrock');

% Barplot of btseq
fig1_bt = figure();
bar(btseq_1);
title(sprintf(t, ft_print));

fig1_surf = figure();
surf(X, Y, Z,'EdgeColor','none');
hold on
plot3([x0(1) xseq_1(1, :)], [x0(2) xseq_1(2, :)], [ros_func(x0), ros_func(xseq_1)], 'r--*');
hold off
title(sprintf(t, ft_print));



