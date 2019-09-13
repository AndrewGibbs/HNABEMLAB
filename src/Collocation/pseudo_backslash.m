function [x,K] = pseudo_backslash(A, b,threshold)
    %function [x,K] = pseudo_backslash(A, B)
    %
    %   Simulate the result of A\B for an underdetermined matrix A.
    %   Implementation based on the description found here:
    %   http://www.mathworks.com/matlabcentral/newsreader/view_thread/238504

    % Compute a reduced QR-decomposition with pivoting
    %[Q,R,E] = qr(A,0);
    [U,S,V] = svd(A,0);

    K = length(find(abs(S)/S(1,1) > threshold));

    S_trunc_inv = zeros(size(S)).';
    S_trunc_inv(1:K,1:K) = (S(1:K,1:K))^(-1);
    x = V*S_trunc_inv*U'*b;
end