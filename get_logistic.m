function [mp, prob, r, t0] = get_logistic( X1, X2, bins )
%GET_LOGISTIC Summary of this function goes here
%   Detailed explanation goes here

%     X1 = P1(tumor_mask);
%     X2 = P1(brain_mask);
%     
v = linspace(0,1,bins + 2);
lb = v(2:bins + 1)';
ub = v(3:bins + 2)';
mp = median([lb ub]')';

for i = 1:bins
    
    prob(i) = sum( (X1 >= lb(i)) & (X1 < ub(i)) ) / ...
        (sum( (X1 >= lb(i)) & (X1 < ub(i)) ) + sum( (X2 >= lb(i)) & (X2 < ub(i)) ));

end

paramguess=[100,0.15];
params=fminsearch(@penalty_logistic,paramguess,[],mp(~isnan(prob)),prob(~isnan(prob)));

r = params(1);
t0 = params(2);

end

