function [ x,z,r, w] = Duffy_quad_genG( a, b, N, width_A)
%See duffy_r.jpg for an explanation of these variables.
%The main aim here is to reduce the rounding error incurred from r=|s-t|,
%by treating the kernel as a function of a single variable r.
  
    if nargin==3
        width_A=b-a; %let the width_A imput be optional.
    end
    
    %compute weights and nodes for L[f]
    [x,w] = quad_gengauss_log(N);
    [gamma,mu] = rearrange(x,x);
    [w1,w2] = rearrange(w,w);
    w_L=w1.*w2;
    %weights and nodes over R[f] will start identically
    w_R=w_L;	sigma=gamma;   zeta=mu;

    W_L=width_A^2.*mu.*w_L;
    x_L=a+mu*width_A;
    z_L=a+(1-gamma).*mu*width_A;
    r_L=gamma.*mu*width_A;

    W_R=w_R.*zeta*width_A^2;
    x_R=a+(1-zeta)*width_A;
    z_R=width_A*(sigma.*zeta+1-zeta)+a;
    r_R=width_A*zeta.*sigma;

    w=[W_L.' W_R.'].';
    r=[r_L.' r_R.'].';
    x=[x_L.' x_R.'].';
    z=[z_L.' z_R.'].';
end
%-------------------------------------------------------------------------%
%Testing in maple against:
%int(int(HankelH1(0, abs(x-y))*x^2*y^3, x = 0 .. 1), y = 0 .. 1.)
%   =0.81962874537748873567e-1-.11487432470092055116i

%error      N
%O(e-16)	10 
%O(e-16)    50 
%0          100