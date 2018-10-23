function x = findNonOscBitR(g,a,b,kwave,C)
%For non-vanishing g', find x such that kwave*|g(x)-g(a)|>=C
    x = a;
    while kwave*abs(g(b)-g(x)) >= C
        x = b - abs(b-x)/2;
    end
end