function [x,w] = gradSingQuad(a, b, N, oscs, layers, delta, abWidth)
%1D brute force quad rule designed to handle singularities and oscillations
    if nargin == 6
        abWidth = b-a;
    end
    [X,W]=gauleg(N);
    X=(X+1)/2; W=W/2;
    wavelengthApprox = abWidth/max(oscs,1);
    %set the initial bit:
    %w = delta^p_max*W*abWidth; x = a+ delta^p_max*X*abWidth;
    w = []; x = [];
    n = layers ;
    %for n = fliplr(1:(layers))
        %if next graded element will be less than half a wavelength...
    while meshElWidth(n) < wavelengthApprox/2 && layers > 0 && n > 0
        %keep grading
        x  = [x; meshElStart(n) + meshElWidth(n).*X];
        w = [w; W*meshElWidth(n)];
        %next layer:
        n = n - 1;
    end
    %end
    
   %fill remainder of interval up with osc resolving meshwidths
   if meshElStart(n) < b
       remainingOscs = ceil((b-meshElStart(n))/wavelengthApprox);
       subWavelengthApprox = (b-meshElStart(n))/remainingOscs;
       for m = 1:remainingOscs
           x = [x; meshElStart(n) + (m-1)*subWavelengthApprox + X*subWavelengthApprox];
           w = [w; W*subWavelengthApprox];
       end
   end
    
    function elWidth = meshElWidth(n)
        if n>layers
            error('requested length of mesh element with index which is too large, does not exist');
        end
        if n == layers
            elWidth = abWidth*delta^(layers-1);
        else
           elWidth = abWidth*delta^(n-1)*(1-delta); 
        end
    end
    
    function elStart = meshElStart(n)
        if n>layers
            error('requested start point of mesh element with index which is too large, does not exist');
        end
        if n == layers
            elStart = a;
        else
            elStart = a + abWidth*delta^(n);
        end
    end
end

