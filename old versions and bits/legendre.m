function P = legendre( n, x )
%   not defined recursively until n=9, as it runs around 35% faster this
%   way.
    switch n
    case 0
        P=ones(size(x));
    case 1
        P=x;
    case 2
        P=1/2*(3*x.^2-1);
    case 3
        P=1/2*(5*x.^3-3*x);
    case 4
        P=1/8*(35*x.^4-30*x.^2+3);
    case 5
        P=1/8*(63*x.^5-70*x.^3+15*x);
    case 6
        P=1/16*(231*x.^6-315*x.^4+105*x.^2-5);
    case 7
        P=1/16*(429*x.^7-693*x.^5+315*x.^3-35*x);
    case 8
        P=1/128*(6435*x.^8-12012*x.^6+6930*x.^4-1260*x.^2+35);
    otherwise %compute using recursive def'n
        P=( (2*n-1)*x.*legendre(n-1,x)-(n-1)*legendre(n-2,x)  )/(n);
    end
end