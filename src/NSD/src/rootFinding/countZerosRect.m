function argPri = countZerosRect( f,df, rect, intThresh, N, moments )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
    if nargin<=4
        N=15;
    end
    
    if nargin<=5
        moments=1;
    end
    
    if length(rect)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    rect(5)=rect(1); %add 5th value as first to complete loop
    
    argPri(1)=.5;
    while (abs(mod(real(argPri(1)),1)-1) > intThresh   &&   abs(mod(real(argPri(1)),1)) > intThresh) || abs(imag(argPri(1))) > intThresh
        %stop while loop getting out of hand:
        if N>10000
            error('Npts to resolve integrand on rectangle edge has become stupidly big');
        end
        %reset output variable
        argPri=zeros(moments,1);
        for j=1:4
                width=abs(rect(j+1)-rect(j));
                quadPts=max(ceil(N*width),N);
                [z_, w] = quad_gauss(quadPts, 0, width);
                %want unit direction of contour along edge
                dir=round((rect(j+1)-rect(j))/width);
                z = rect(j) + dir*z_;

            for m=1:moments
                argPri(m)=argPri(m)+w.'*dir*(z.^(m-1).*df(z)./f(z))/(2*pi*1i);
            end
        end
        N=N*2;
    end
    
end