function isZero = isZeroOnLine( a,b,F, dF, bigThresh, lineLength )

    if nargin<=5
        lineLength=abs(b-a);
    end

    lineMap=@(s) a+exp(angle(b-a)*1i)*lineLength.*s;
    
    f=@(s) F(lineMap(s));
    df=@(s) dF(lineMap(s));
    
    t=linspace(0,1,15*ceil(max(1,lineLength)));
    
    %the following will blow up near a pole
    if max(abs(df(t)./f(t)))>bigThresh
        isZero=true;
    else
        isZero=false;
    end
end