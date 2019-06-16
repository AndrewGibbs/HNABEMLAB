function [ G, SPs, SPOs ] = SqrtTypePhase( a,b )
%for phase functions of type g(x)=sqrt(x^2+a)+bx

%returns each derivative of g(x) (as many as will be required by NSD code)
%alongside stationary points and orders

%could ultimately change this to output inverses as well... if they're not
%too tricky to find...

zeroThresh=1E-15;

    %deal with this special case first;
    if abs(a)<zeroThresh
        G={@(x) (1+b)*x, @(x) (1+b), @(x) 0};
        SPs=[]; SPOs=[];
        return;
    end
    
    %excliding the above special case, can define the derivatives more generally:
    G={@(x) (x.^2+a).^.5 + b*x, @(x) x./sqrt(x.^2+a) + b, @(x) a./(a+x.^2).^(3/2)};
    
    %and a second special case:
    if abs(b)<zeroThresh
        SPs=0; SPOs=1;
    elseif abs(b^2-1)<zeroThresh
        SPs=[]; SPOs=[];
    else
        %test both possible stationary points, they don't always work:
        SPs=[]; SPOs=[];
        if abs(G{2}(sqrt(a*b^2/(1-b^2)))) < zeroThresh
            SPs=[SPs sqrt(a*b^2/(1-b^2))];
        end
        if abs(G{2}(-sqrt(a*b^2/(1-b^2)))) < zeroThresh
            SPs=[SPs -sqrt(a*b^2/(1-b^2))];
        end
        SPOs=ones(size(SPs));
    end
    
end


