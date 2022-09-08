function [gradfx] = findiff_grad(f, x, h, type)
%
% function [gradf] = findiff_grad(f, x, h, type)
%
% Function that approximate the gradient of f in x (column vector) with the
% finite difference (forward/centered) method.
%
% INPUTS:
% f = function handle that describes a function R^n->R;
% x = n-dimensional column vector;
% h = the h used for the finite difference computation of gradf
% type = 'fw' or 'c' for choosing the forward/centered finite difference
% computation of the gradient.
%
% OUTPUTS:
% gradfx = column vector (same size of x) corresponding to the approximation
% of the gradient of f in x.

gradfx = zeros(size(x));

switch type
    case 'fw'
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x))/ h;
            % ALTERNATIVELY (no xh)
            % gradf(i) = (f([x(1:i-1); x(i)+h; x(i+1:end)]) - f(x))/h;
        end
    case 'c'
        for i=1:length(x)
            xh_plus = x;
            xh_minus = x;
            xh_plus(i) = xh_plus(i) + h;
            xh_minus(i) = xh_minus(i) - h;
            gradfx(i) = (f(xh_plus) - f(xh_minus))/(2 * h);
            % ALTERNATIVELY (no xh_plus and xh_minus)
            % gradf(i) = ((f([x(1:i-1); x(i)+h; x(i+1:end)]) - ...
            %     f([x(1:i-1); x(i)-h; x(i+1:end)]))/(2*h);
        end
    otherwise % repeat the 'fw' case
        for i=1:length(x)
            xh = x;
            xh(i) = xh(i) + h;
            gradfx(i) = (f(xh) - f(x))/h;
            % ALTERNATIVELY (no xh)
            % gradf(i) = (f([x(1:i-1); x(i)+h; x(i+1:end)]) - f(x))/h;
        end
end



end

