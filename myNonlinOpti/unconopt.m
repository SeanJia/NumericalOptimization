function [x,fval] = unconopt(f, x0)

% This function use modified Newton method with backtracking to compute
% the real vector x for given function func, where f(x) achieves 
% the local minimizer around the starting point x0;
%
% Usage:
% input: f: the given function; this should be declared as the example
%           below in command line, for f(x,y,z) = x^2+y^2+z^2, type:
%                         f = inline('x^2+y^2+z^2');
%        x0: the initial guess, with default value as the origin
% output: x: the resulting x
%        fval: the resulting f(x)
%
% Author: Zhiwei Jia

args = argnames(f);                  % the variables
x = zeros(length(args),1);           % the initial guess
tol = 1.0e-10;                       % tolerence for zero finding
fx = sym(char(f));                   % convert to symbolic version
gx = jacobian(fx)';                  % compute its gradient
Hx = jacobian(gx);                   % compute its Hessian

% the valid inital guess
if length(x0) ~= length(args)
    x0 = x;
end

% strings for f,g, and H
f = char(fx);
g = char(gx);
H = char(Hx);

fprintf('\n  Unconstrained optimization for the function %s \n', f);
[x,fval] = newton(f,g,H,x0,tol);   % compute the value
fprintf('\n  The final result:\n  f is %f, and x is\n', fval);
fprintf('                           %f\n',x);

% if more than one variable, display their name in order
if length(args) > 1
    fprintf('  where the vector x stands for the variable\n');
    for i = 1:length(args)
        disp(args(i));
    end
end
end