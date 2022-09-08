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
h = sqrt(eps);
% ROSENBROCK FUNCTION
ros_func = @(x)100*(x(2,:)-x(1,:).^2).^2+(1-x(1,:)).^2;
ros_grad = @(f,x,h,type) findiff_grad(f, x, h, type);
ros_hess = @(f, x, h)findiff_Hess(f, x, h);

x0 = [1.2; 1.2];
%x0 = [-1.2; 1];
%x0 = [-5; -5];
disp("x0: " + mat2str(x0));

disp("****** TESTING ROSENBROCK FUNCTION *******");
disp("**** NEWTON METHOD WITH 'FW' FINITE DIFERENCE  *****");

tic;  
[xk_1, fk_1, gradfk_norm, k_1, xseq_1, btseq_1] = ...
    newton_general(x0, ros_func, ros_grad, ros_hess, kmax, ...
    tolgrad, c, rho, btmax, 'fw', 'c', h);
ex = toc;
disp('************** FINISHED ***************');
disp('************** RESULTS ****************');
disp(['xk: ', mat2str(xk_1), ' (actual minimum: [1; 1]);']);
disp(['f(xk): ', num2str(fk_1), ' (actual min. value: 0);']);
disp(['N. of Iterations: ', num2str(k_1),'/',num2str(kmax), ';']);
fprintf('EXCUTION TIME IS : %.2f ms \n' , ex *100 );

disp('************************************');

disp("---");


% PLOTS
% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 25, 500));
Z = 100*(Y-X.^2).^2+(1-X).^2;

fig1_n = figure();
% Contour plot with curve levels for each point in xseq
[C, ~] = contour(X, Y, Z);
hold on;
% plot of the points in xseq
plot([x0(1) xseq_1(1, :)], [x0(2) xseq_1(2, :)], '--*');
hold off;
title("NEWTON METHOD WITH 'FW' FINITE DIFERENCE ' - Rosenbrock" );

% Barplot of btseq
fig1_bt = figure();
bar(btseq_1);
title("NEWTON METHOD WITH 'FW' FINITE DIFERENCE ' - Rosenbrock");

fig1_surf = figure();
surf(X, Y, Z,'EdgeColor','none')
hold on
plot3([x0(1) xseq_1(1, :)], [x0(2) xseq_1(2, :)], [ros_func(x0), ros_func(xseq_1)], 'r--*')
hold off
title("NEWTON METHOD WITH 'FW' FINITE DIFERENCE ' - Rosenbrock")




