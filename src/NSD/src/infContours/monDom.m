function R = monDom(a)
%returns an upper bound on the radius at which the leading term of a
%mononial becomes the dominant term
%a should be an (N+1)x1 vector of coefficients such that
%p(x) = a_0 + a_1x + a_2x^2+... a_Nx^N
    if a(end)~=0
        aNormal = a(1:(end-1))/a(end) ;
    else
        aNormal = a ;
    end
    aMax = max(1,max(aNormal(2:end)));
    R = aMax + a(1);
end

