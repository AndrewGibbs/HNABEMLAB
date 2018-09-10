function x = findNonOscBit(g,a,b,kwave,C)
%For non-vanishing g', find x such that kwave*|g(x)-g(a)|>=C
    x = b;
    while kwave*abs(g(x)-g(a)) >= C
        x = a + abs(x-a)/2;
    end
end

