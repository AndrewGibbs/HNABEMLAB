function [x,w] = LaplaceSD(start, freq, n, singular)
%Ureturns weights and nodes for evaluating Laplace-type SD integrals in
%complex plane
    if singular
        [x_,w_] = quad_gengauss_loglaguerre(2*n);
    else
        [x_, w_] = gausLagHC(n); %use hardcoded values for speed
        %[x_,w_] = quad_gauss_exp(1, n);
    end
    x = start + 1i*x_/freq;
    w = 1i*exp(1i*freq*start)*w_/freq;
end