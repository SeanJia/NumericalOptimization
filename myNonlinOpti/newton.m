function [x,fval] = newton(f,g,H,x0,tol)

% This function use newton's method to find a local minimizer of a nonlinear
% function f(x).
%
% Usage:
% input: f: a string containing the function's expressiong
%        g: a string for the gradient of f
%        H: a string for the Hession
%        x0: the initial guess tarting guess
%        tol: the tolerance for considering zero
% output: x: the final x
%         fval: the final result of the function
%
% Author: Zhiwei Jia

itmax = 50;             % limit on function evaluations
iter = 0;               % number of iterations           
x = x0;                 % the initial guess
dpmax = 20;             % limit on damping steps
mu = 1/4;               % sufficient decrease paramter
delta = 1.0e-10;        % smallest eigenvalue cut-off
bt = 0;                 % total backtracking iterations

args = argnames(inline(f));  % the variables

% f,g and H in symbolic version
f = sym(f);                 
g = sym(g);
H = sym(H);

% evaluation of f,g and H
fval = double(subs(f,args,x));
gval = double(subs(g,args,x));
Hval = double(subs(H,args,x));

% when initial guess has gradient zero, try surrounding points
while double(gval) < tol
    x = x + 0.1;
    gval = double(subs(g,args,x));
end

% print statistics
fprintf('  Itn   bt   Step        f(x)       norm(g)\n');
str1 = sprintf(' %3g %4g %9.2e', iter, bt, 0);
str2 = sprintf(' %14.7e %9.2e' , double(fval), norm(double(gval)));
disp([str1 str2]);  

% the while loop for newton iteration
while iter < itmax && norm(gval) > tol 
    
    % modified eigenvalue
    [V,D] = eig(Hval);
    for i = 1:length(D)
        D(i,i) = max([delta, abs(D(i,i))]);
    end
    
    % obtain the newton direction
    p = -V*(inv(D))*V'*double(gval);
    
    % backtracking part
    alpha = 1;
    it = 0;
    newx = x + alpha*p;
    newf = double(subs(f,args,newx));
    while newf >= fval + alpha*mu*gval'*p && it < dpmax
        alpha = alpha/2;
        it = it + 1;
        newx = x + alpha*p;
        newf = double(subs(f,args,newx));
    end
    
    % update the relevant info
    bt = bt + it;
    x = newx;
    iter = iter + 1;
    fval = double(subs(f,args,x));
    gval = double(subs(g,args,x));
    Hval = double(subs(H,args,x));
    str1  = sprintf(' %3g %4g %9.2e', iter, bt, alpha);
    str2  = sprintf(' %14.7e %9.2e', fval, norm(gval));
    disp([str1 str2])
end
end
