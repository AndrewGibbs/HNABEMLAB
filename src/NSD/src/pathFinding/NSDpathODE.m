function dHdp = NSDpathODE(p, h, n, g, ICs, ascFlag, pSmallThresh, hTaylor)
    %computes NSD path h and h', at point order n

    AscDescDir=1; %
    
    if ascFlag
        AscDescDir=-1;
    end
        
    if nargin <=6
        pSmallThresh=1E-12;
    end
    dgSmallThresh = pSmallThresh;
    
    if nargin<=7
        noTaylorFlag = true;
    else
        noTaylorFlag = false;
    end
    noTaylorFlag = true;
    if length(g)<(2*n)
        error('Higher derivatives of g required as input');
    end
    
    if n==0
         dHdp =  AscDescDir*1i./g{2}(h); 
    elseif p>=pSmallThresh && abs(g{2}(h(1)))>=dgSmallThresh
        switch n
            case 1
                dHdp = [h(2); (AscDescDir*2i-h(2).^2.*(g{3}(h(1))))./g{2}(h(1))]; 
                % H=[h,h']

            case 2
                dHdp = [h(2); h(3); (AscDescDir*6i-h(2).^3.*(g{4}(h(1)))-3*g{3}(h(1)).*h(2).*h(3))./g{2}(h(1))];
                % H=[h,h',h'']
            case 3
                dHdp = [h(2); h(3); h(4); (AscDescDir*24i - g{5}(h(1)).*h(2).^4 - 6*g{4}(h(1)).*h(2).^2.*h(3) - g{3}(h(1)).*(3*h(3).^2 + 4*h(4).*h(2)))./g{2}(h(1))];
                % H=[h,h',h'', h''']
            case 4
                dHdp = [h(2); h(3); h(4); h(5); ...
                        (AscDescDir*120i - g{6}(h(1)).*h(2).^5 - 10*g{5}(h(1)).*h(2).^3.*h(3) - 10*h(4).*g{4}(h(1)).*h(2).^2 - h(3).*g{3}(h(1)) - 5*h(2).*(3*g{4}(h(1))).*h(3).^2 + h(5).*g{3}(h(1)))./g{2}(h(1));
                        ];
                
            otherwise
                error('Cant go above 3rd order stationary points yet');
        end
    else
        %very small p, can get unstable... so approximate Taylor style:
        %dHdp = [ICs(2); zeros(n,1)];
        if noTaylorFlag
            dHdp = [ICs(2:(n+1)).'; 0];
        else
            for n_=1:(n+1)
                dHdp(n_) = hTaylor.path(p,n_);
            end
            dHdp=dHdp(:);
        end
    end
end