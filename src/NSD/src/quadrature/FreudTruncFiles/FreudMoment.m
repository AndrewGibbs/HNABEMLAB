function m = FreudMoment( b, p, n )
%computes j'th moment of Freud-type weight w(x)=exp(-x^p) for p=1,2,... over inverval
%[0,b]
%result follows by substituting s=t^p into the integral of t^n*exp(-t^p)
%over [0,b], and by subsequently acknowledging that Mathworks' definition of the
%incomplete gamma function is STUPID

    a=(1+n)/p;
    m=gamma(a).*gammainc(b^p,a,'lower')/p;

end