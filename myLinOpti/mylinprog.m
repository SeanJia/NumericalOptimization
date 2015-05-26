function [cval,x] = mylinprog(A,b,c)

% This function use Simplex method for linear programming
% in all inequality form; it compute the following  
% linear optimizatin problem:
%            minimize     c'x
%           subject to   Ax >= b,  x >= 0
% Usage: 
%        input: A, b: relevant constraint info
%               c: objective function to be minimized
%        output: x: optimal point
%                objvalue: resulting objective value
%
% During each iteration, following message will be displaieded:
%     Itn: the current iteration number.
%     Objective: the value of the current objective function.
%     MinofLM: the Lagrange multiplier of the constraint selected
%              to leave the working set at the start of the
%              present iteration.
%     Del: the constraint selected to leave the working set
%     Add: the constraint selected to enter the working set
%     Step: the steplength taken along the current search 
%
% Author: Zhiwei Jia

[m,n] = size(A); 
A0 = [A ones(m,1)];       % constraint matrix for phase 1 lin prog
A0(find(b<=0),n+1) = 0;   % add the new variable theta
c0 = [zeros(n,1);1];      % obj function for phase 1 lin prog
b0 = b; 
theta0 = max(max(b),0); 
x0 = [zeros(n,1); theta0];    % the initial vertex for phase 1 LP
fprintf('\n  Phase 1 LP for finding an initial vertex:\n\n');
[cval,x] = simplex(A0,b0,c0,x0);
fprintf('\n  Phase 2 LP for the initial optimization problem:\n\n');
[cval,x] = simplex(A,b,c,x(1:n,:));

end

