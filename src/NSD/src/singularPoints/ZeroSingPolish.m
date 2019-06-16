function [zeros, sings, zeroCount, singCount] = ZeroSingPolish( points,  f, df, thresh, N )
%takes estimated values of zeros/singularities, deletes points which are too close to others,
%then recomputes more acurately, dividing by counting function

zeros=[]; sings=[]; zeroCount=[]; singCount=[];

    %first detect 'repeated' (up to some distance threshold) points
    keepIndices=true(size(points));
    for p=1:length(points)
        for q=(p+1):length(points)
           if abs(points(p)-points(q))<thresh %have found the same point twice
               keepIndices(q)=false;
           end
        end
    end
    
    %now delete repeated* values
    points=points(keepIndices);

    for p=1:length(points)
        %point definately lies inside this square:
        square=points(p) + thresh*[1+1i  -1+1i  -1-1i  1-1i];
        %get counts and counts*zeros
        countOut=countZerosSingsRect( f, df, square, N, 2 );
        count=countOut(1);
        count_X_z0=countOut(2);
        
        %now determine type of point, and compute acurately it's value
        if count>0 %it's a zero
            zeros=[zeros count_X_z0/count];
            zeroCount=[zeroCount count];
        else %it's a singularity
            sings=[sings count_X_z0/count];
            singCount=[singCount count];
        end
        clear square countOut count count_X_z_0;
    end

end
