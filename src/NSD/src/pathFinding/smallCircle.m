function C = smallCircle(c,r,N)
    x = linspace(0,2*pi,N);
    C = c + r*cos(x) + r*1i*sin(x);
end