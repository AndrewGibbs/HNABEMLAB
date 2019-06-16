function [X,W] = infGradedQuad(N, p_max, delta , m)
%approximates infinite Freud-type integral by finite integral, so that can
%we can grade torwads zero, to account for singularities. Truncating
%Freud-type integrals for m>1 is not possible in the same way as m=1.
    thresh = eps;
    L = 50 *(-log(thresh))^(1/m);
    [ x, w ] = GradedQuad( N, p_max, delta );
    X = L*x;
    W = L*w.*exp(-(X.^m));
end

