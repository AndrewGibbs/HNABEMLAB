function [x,K] = pseudo_backslash(A, B,threshold)
%function [x,K] = pseudo_backslash(A, B)
%
%   Simulate the result of A\B for an underdetermined matrix A.
%   Implementation based on the description found here:
%   http://www.mathworks.com/matlabcentral/newsreader/view_thread/238504

% Compute a reduced QR-decomposition with pivoting
[Q,R,E] = qr(A,0);

% Find the numerical rank K by truncating R
%threshold = 2e-13;  % value determined experimentally
Rd = diag(R);
L = length(find(abs(Rd) < threshold));
if L > 0
    K = find(abs(Rd) < threshold, 1 ) - 1;
else
    K = length(Rd);
end

K;

% Select the first K columns on R and solve for K unknowns
R1 = R(:,1:K);
x1 = R1 \ (Q'*B);

% Copy the values of x1 to their right place in the solution vector, using
% the permutation returned by the pivoted QR
x = zeros(size(A,2),1);
x(E(1:K)) = x1;
