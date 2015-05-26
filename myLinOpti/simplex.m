function [cval,x] = simplex(A,b,c,x0)

% Simplex method for linear optimizatin problem:
%     minimize       c'x
%     subject to   Ax >= b,   x >= 0
%     where x0 is the initial vertex.
%
% Author: Zhiwei Jia

unbounded = 0;   
maxStep = 1;
tol = 1.0e-8;   % the tolerence for consider a number as zero
p = 1;          % initial search direction
x = x0;         % initial vertex
itnMax = 50;    % maximum iteration number
itn = 0;        % current iteration
[m,n] = size(A);          % size of the constraint matrix
Afull = [A; eye(n)];      % full version of constraint matrix 
bfull = [b; zeros(n,1)];  % full version of vector b
r = Afull*x - bfull;      % the residual
vio_const = r < -tol*ones(m+n,1);   % indexes of violated constraints 
status_const = ones(m+n,1);         % status of each constraint
active_const = find(abs(r) < tol);  % indexes of active constraints
status_const(active_const) = 0;     
Aw = Afull(active_const,:);   % the working set of constraints
bw = bfull(active_const);     % corresponding part for b
cval = 0;                     % the objective funcion value

% when the given initial point is not feasible
if any(vio_const) 
    error('ERROR! Initial point is infeasible!');
end

% when the given initial point is not a vertex
C = active_const'; 
if length(active_const) > n
    C = nchoosek(active_const,n); 
end
[many,n] = size(C);  % the number of possible combination
notvertex = 1;       % whether the point is a vertex
for i = 1:many       % check for each combination
    active_const = C(i,:);
    Aw = Afull(active_const,:); 
    if rank(Aw) == n
        notvertex = 0;
    end
end
if notvertex
    error('ERROR! The initial point is not a vertex!');
end

fprintf('  Simplex method for min c''*x, subject to A*x >= b, x >=0\n');
fprintf('  Itn      Objective     minofLM   Del  Add      Step\n');

% while loop for the simplex iteration
while itn < itnMax && norm(maxStep*p) > tol
    
    % find the Lagrange multiplier and its most nagative component
    % to check for optimality
    lambda = Aw'\c;   
    [minval,minidx] = min(lambda);
    fprintf(' %3d   %14.7f  %6d', itn, c'*x, minidx);
    
    % the optimal case
    if minval >= 0
        fprintf('    Optimal\n');
        break; 
    end

    % a standard vector whose minidx position is 1
    e_minidx = zeros(n,1); 
    e_minidx(minidx) = 1;
    
    % the descent direction
    p = Aw\e_minidx;
    fprintf('    %4d ', active_const(minidx)); % the constraint to leave
                                               % the working set

    % obtain the max step and relevant info for the desccent direction                                           
    inactive = find(r > 0);
    [maxStep,i,unbounded] = maxstep(r,Afull*p,inactive);
    % the unbounded case
    if unbounded
        fprintf('  Unbounded step\nso this problem is unbounded below');
        unbounded = 1;
        break; 
    end     
    
    status_const(i) = 0;                      % becoming active                    
    status_const(active_const(minidx)) = 1;   % becoming inactive
    active_const = find(status_const == 0);   % update active set
    Aw = Afull(active_const,:);               % update working matrix
    bw = bfull(active_const,:);               % corresponding vector b
    x = x + maxStep*p;                            % the new vertex
    r = Afull*x - bfull;                      % new residual
    fprintf('%4g      %8.2e\n', i, maxStep);
    itn = itn + 1;
    cval = c'*x;
end

% show the result
if ~unbounded
    fprintf('\n  The final objective is %f, the final x is \n',cval);
    fprintf('      %f\n', x);
end
end
