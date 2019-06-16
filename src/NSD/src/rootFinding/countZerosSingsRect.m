function argPri = countZerosSingsRect( f,df, rect, N, moments )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
    if nargin<=3
        N=15;
    end
    
    if nargin<=4
        moments=1;
    end
    
    if length(rect)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    rect(5)=rect(1); %add 5th value as first to complete loop
    
    argPri=zeros(moments,1);
    for j=1:4
            width=abs(rect(j+1)-rect(j));
            quadPts=max(ceil(N*width),5);
            [z_, w] = quad_gauss(quadPts, 0, width);
            %want unit direction of contour along edge
            dir=round((rect(j+1)-rect(j))/width);
            z = rect(j) + dir*z_;

        for m=1:moments
            argPri(m)=argPri(m)+w.'*dir*(z.^(m-1).*df(z)./f(z))/(2*pi*1i);
        end
    end
    
end