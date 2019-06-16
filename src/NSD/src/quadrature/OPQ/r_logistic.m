% R_LOGISTIC Recurrence coefficients for the logistic weight
% function.
%
%    ab=R_LOGISTIC(n) generates the first n recurrence coefficients
%    for the weight function w(t)=exp(-t)/(1+exp(-t))^2.
%
function ab=r_logistic(N)
if N<=0, error('N out of range'), end
mu=1; n=(1:N-1)';
ab=[zeros(N,1) [mu;n.^4*pi^2./(4*n.^2-1)]];
