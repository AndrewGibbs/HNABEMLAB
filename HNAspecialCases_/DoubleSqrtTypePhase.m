function [ G, SPs, SPOs ] = DoubleSqrtTypePhase( a, b, c, d )
%for phase functions of type g(x)=sqrt(x^2+ax+b)+sqrt(x^2+cx+d)

    zeroThresh=1E-15;

    G{1}=@(x) sqrt(x.^2+ax+b)+sqrt(x.^2+cx+d);
    G{2}=@(x) 1/2*((a + 2*x)./sqrt(x.*(a + x) + b) + (c + 2*x)./sqrt(x.*(c + x) + d));
    G{3}=@(x) 1/2*(-(a + 2*x).^2./(2*(x*(a + x) + b).^(3/2)) + 2./sqrt(x*(a + x) + b) + 2./sqrt(x*(c + x) + d) - (c + 2*x).^2./(2*(x*(c + x) + d).^(3/2)));
    
    
    %SPs=roots(d+(a/2)^2-b-(c/2)^2 , a*d+c-b*c-a, d*(a/2)^2-b*(c/2)^2);
    % symbolic me, Tuesday

    %SPs=roots([a^2-c^2+4*a^2*c-4*a*c^2+4*d-4*b,   a^2*c-c^2*a+4*a*d-4*b*c,   a^2*d+c^2*b]);
    %symbolic me, Wednesday

    %SPs=-(b - d)/(a - c); % symbolic matlab - wrong

    % wolfram alpha: first root seems wrong, second OK. But why...
    SPs(1)=(-sqrt(a^2 - 4*b)*(a - c)*sqrt(c^2 - 4*d) + a^2*c - a*c^2 + 4*a*d - 4*b*c)/(2*(-a^2 + 4*b + c^2 - 4*d));

    SPs(2)=(sqrt(a^2 - 4*b)*(a - c)*sqrt(c^2 - 4*d) + a^2*c - a*c^2 + 4*a*d - 4*b*c)/(2*(-a^2 + 4*b + c^2 - 4*d));
    
    copySPs=SPs;
    SPs=[]; SPOs=[];
    for j=[1 2]
       if abs(G{2}(copySPs(j)) )<zeroThresh
           SPs=[SPs copySPs(j)];
           if abs(G{3}(copySPs(j)) )<zeroThresh
               SPOs=[SPOs 2];
           else
               SPOs=[SPOs 1];
           end
       end
    end
    
end
