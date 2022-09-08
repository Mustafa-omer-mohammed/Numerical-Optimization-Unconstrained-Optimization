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

% ROSENBROCK FUNCTION
ros_func = @(x)100*(x(2,:)-x(1,:).^2).^2+(1-x(1,:)).^2;
ros_grad = @(x)[400*x(1,:).^3-400*x(1,:).*x(2,:)+2*x(1,:)-2; ...
    200*(x(2,:)-x(1,:).^2)];
ros_hess = @(x)[1200*x(1)^2-400*x(2)+2,-400*x(1);-400*x(1),200];

x0 = [1.2; 1.2];
% x0 = [-1.2; 1];
%x0 = [-5; -5];
disp("x0: " + mat2str(x0));



disp("*** STEEPEST DESCENT WITH BACKTRACKING **");
tic
[xk_2, fk_2, gradfk_norm_2, k_2, xseq_2, btseq_2] = steepest_descent_bcktrck(x0, ros_func, ...
    ros_grad, alpha0, 1000, 1e-8, c, rho, btmax,1);
ex = toc;
disp('************** FINISHED ***************');
disp('************** RESULTS ****************');
disp(['xk: ', mat2str(xk_2), ' (actual minimum: [1; 1]);']);
disp(['f(xk): ', num2str(fk_2(end)), ' (actual min. value: 0);']);
disp(['N. of Iterations: ', num2str(k_2),'/',num2str(10000), ';']);
fprintf('EXCUTION TIME IS : %.2f ms \n' , ex *100 );
disp('************************************');

% line plot of function value and gradient norm
figure();
yyaxis left;
semilogy(fk_2, 'LineWidth', 2);
ylabel('Value of the function');
yyaxis right;
semilogy(gradfk_norm_2, '--', 'LineWidth', 2);
ylabel('Norm of the gradient');
legend('Fk trend', 'GradFk trend', 'Location', 'southwest');
xlabel('Number of iteration (k)');
disp('---');
% PLOTS
% Creation of the meshgrid for the contour-plot
[X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 25, 500));
Z = 100*(Y-X.^2).^2+(1-X).^2;


fig2_n = figure();
% Contour plot with curve levels for each point in xseq
[C, ~] = contour(X, Y, Z);
hold on;
% plot of the points in xseq
plot([x0(1) xseq_2(1, :)], [x0(2) xseq_2(2, :)], '--*');
hold off;
title('Steepest descent with backtracking - Rosenbrock');

% Barplot of btseq
fig2_bt = figure();
bar(btseq_2);
title('Steepest descent with backtracking - Rosenbrock');

fig2_surf = figure();
surf(X, Y, Z,'EdgeColor','none')
hold on
plot3([x0(1) xseq_2(1, :)], [x0(2) xseq_2(2, :)], [ros_func(x0), ros_func(xseq_2)], 'r--*')
hold off
title('Steepest descent with backtracking - Rosenbrock')


