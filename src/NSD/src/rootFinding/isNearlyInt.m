function yn = isNearlyInt( x,T )
%determines if x is nearly an integer, up to some threshold T
    xReal=real(x);
    xImag=imag(x);
    
    if abs(xImag)>T
        yn=false;
        return;
    end
    
    r=mod(xReal,1); %remainder when real part is divided by one
    if r>=.5
        r=1-r;
    end
    
    if r<=T
        yn=true;
    else
        yn=false;
    end

end

