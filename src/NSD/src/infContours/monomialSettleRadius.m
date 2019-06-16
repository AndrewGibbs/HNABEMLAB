function R = monomialSettleRadius(a, Npts, freq, thresh)
%estimates the radius of a ball centred at the origin (perhaps
%unecessarily), outside of which a monomial is dominated by it's leading
%term, such that the steepest descent paths follow the valleys closely

    if nargin<=2
        Npts=15;
    end
    if nargin<=3
        freq=1; %extra freq info can only speed up process
    end
    if nargin<=3
        thresh=1E-12; % cost of getting from circle edge to infinity
    end
    reScale=sqrt(2);  % amount to rescale by if the circle is too small
    
    
    %reconstruct the phase:
    g = @(z) polyval(a,z);
    
    %get initial guess at radius:
    R = monDom(a);
    
    while thresh < inTheValleys(a, g, R, Npts, freq) %whilst it fails the circle test
        R=R*reScale;
    end
    
end

