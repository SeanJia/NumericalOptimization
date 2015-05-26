function [x,fval] = myzerofinder(F)

% This function use modified Newton method with backtracking to compute
% the real vector x for given function func, where F(x) = 0;
%
% Usage:
% input: F: the given function; this should be declared as the example
%           below in command line, for F(x,y,z) = x^2+y^2+z^2, type:
%                         F = inline('x^2+y^2+z^2');
% output: x: the resulting x
%        fval: the resulting F(x)
%
% Author: Zhiwei Jia

args = argnames(F);                  % the variables
x0 = zeros(length(args),1);          % the initial guess
tol = 1.0e-10;                       % tolerence for zero finding
fx = 0.5*(sym(formula(F)))^2;        % f = 1/2||F(x)||^2, symbolic version       
gx = jacobian(fx)';                  % compute its gradient
Hx = jacobian(gx);                   % compute its Hessian

% strings for f,g, and H
f = char(fx);                        
g = char(gx);
H = char(Hx);

fprintf('\n  Finding zero of the function %s \n', formula(F));
fprintf('  This problem is equivalent to the unconstrained optimization for %s \n\n',char(fx));
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


                     