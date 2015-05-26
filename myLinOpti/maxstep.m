function [maxStep,i,unbounded] = maxstep(r,Ap,inact)

% Given a feasible point  x  for the linear inequalities
% Ax >= b,  maxstep computes maxStep, the maximum step to the
% nearest constraint along the given direction p.
%
% Input: r: the residual Ax - b
%        Ap: the vector A*p
%        inactive: indixes of the inactive constraints
%
% Output: maxStep: the maximum step to the nearest decreasing constraint
%         i: index of such constraint; if more than one blocking constraint
%            return the one with the most negative component of Ap
%         unbounded: whether the solution is unbounded below
%
% Author: Zhiwei Jia

tol = 1.0e-8;         % tolerence of a value to be considered as zeroe
inf = 1.0e10;         % considred as the infinite step
Apinact = Ap(inact);  % the inactive component of A*p
rinact = r(inact);    % the inactive component of residual
decre = find(Apinact < 0);  % indexes of the decreasing constraints

% when there is no dcreasing constraint, the solution is unbounde below
if isempty(decre)
    unbounded = 1;
    maxStep = inf;
    i = 0;
else
    unbounded = 0;
    stepRatios = rinact(decre)./(-Apinact(decre));
    [maxStep,~] = min(stepRatios);
    blocking = find(stepRatios <= maxStep + tol);  % indexes of blocking constraints
    [~,i2] = max(-Apinact(decre(blocking)));       % the one with most nagative Ap
    i = blocking(i2);
    maxStep = stepRatios(i);
    i = inact(decre(i));
end
end
