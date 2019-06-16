function [zeros, sings, zeroCount, singCount] = findZerosSingsRect( f, df, rectIn, thresh, N, info )
%detects roots of analytic function inside of a given rectangle in complex
%domain
%N is number of points per unit length of rectangle
%thresh is the width of the rectangle which we pinpoint the zeros to

MaxMoment=4; %take this many moments as safety net for if the order of a zero cancels with the order of a pole
zeros=[]; sings=[]; zeroCount=[]; singCount=[];
countThresh=.01; %order of singularities/stationary points will be ignored beyond this
bodgeCount=0;

    if nargin<=4
        N=15;
    end
    
    if length(rectIn)~=4
        error('rectangle must consist of four complex values, ordered anti-clockwise');
    end
    
    if nargin<=5
        info=false;
    end
    
    height=max(imag(rectIn))-min(imag(rectIn));
    width=max(real(rectIn))-min(real(rectIn));
    
    
    %first check if zeros lie on edge of rectangle
%     fullRectIn=[rectIn rectIn(1)];
%     for j=1:4
%         zeroOnLine=~isempty(isZeroOnLine( fullRectIn(j), fullRectIn(j+1), f ));
%         if zeroOnLine %if they lie on edge of rectangle
%             %retry with a larger rectangle
%             newRect = rectIn + (rectIn - mean(rectIn))/10 ;
%             allZerosAndSings = findZerosRect( f, df, newRect, thresh, N );
%             return;
%         end
%     end
    
    %first guess, to kick off the while loop
    allZerosAndSings=rand+rand*1i;
    rect=rectIn;
    
    while max(height,width)>thresh % max(abs(f(allZeros)))>thresh %was an alternative
        
        BLcorner=min(real(rect))+1i*min(imag(rect));

        bisectRatio=.5; %start by bisecting in the middle
        
        %zeroOnLine=true; %this value keeps choosing new bisection point until bisection line doesn't land on a zero
        %while zeroOnLine
            if width>height
                %take bototm left corner, and add half (width) rectangle to it
                newRectA=BLcorner + [0  bisectRatio*width  (bisectRatio*width+1i*height)  1i*height];
                newRectB=newRectA + bisectRatio*width;
%                 bisectingLine=BLcorner+bisectRatio*width+[0  1i*height];
%                 bisectingLineLength=height;
            else
                %take bototm left corner, and add half (height) rectangle to it
                newRectA=BLcorner + [0  width  (width+bisectRatio*1i*height)  bisectRatio*1i*height];
                newRectB=newRectA + bisectRatio*1i*height;
%                 bisectingLine=BLcorner+bisectRatio*1i*height+[0  width];
%                 bisectingLineLength=width;
            end
            
%             zeroOnLine=~isempty(isZeroOnLine( bisectingLine(1), bisectingLine(2), f, bisectingLineLength ));
%             bisectRatio=.5+.5*(rand-.5);% in (.25,.75)
%         end
        if info
            plotRects( newRectA, newRectB );
        end
        %count the number of zeros/singularities in each subdivision:
        Acount=max(abs(countZerosRect( f, df, newRectA, N, MaxMoment )));
        Bcount=max(abs(countZerosRect( f, df, newRectB, N, MaxMoment )));
        
        if (Acount)>countThresh && abs(Bcount)>countThresh
            %call self recursively over each sub-tangle, until zeros are
%             %located
%             [zerosA, singsA, zeroCountA, singCountA]=findZerosSingsRect( f, df, newRectA, thresh, N );
%             [zerosB, singsB, zeroCountB, singCountB]=findZerosSingsRect( f, df, newRectB, thresh, N );
            [zerosA, singsA]=findZerosSingsRect( f, df, newRectA, thresh, N, info );
            [zerosB, singsB]=findZerosSingsRect( f, df, newRectB, thresh, N, info );
            zeros=[zerosA zerosB]; sings=[singsA singsB];
            % zeroCount=[zeroCountA zeroCountB]; singCount=[singCountA singCountB];
            %polish things up with this function:
            [zeros, sings, zeroCount, singCount] = ZeroSingPolish([zeros sings],  f, df, thresh, N^2 );
            return;
        elseif (Acount)>countThresh
            rect=newRectA;
%             allZerosAndSings=mean(rect);
        elseif (Bcount)>countThresh
            rect=newRectB;
%             allZerosAndSings=mean(rect);
        elseif (Acount)<=countThresh && (Bcount)<=countThresh
%             allZerosAndSings=[];
%             break
            return;
        else
            if bodgeCount>10
                error('Contour integration has returned a non-complex value over ten times');
            end
            bodgeCount=bodgeCount+1;
            %try shifting the rectangle by a random amount
            rect=mean(rect)+1.1*normaliseRectangle( rect, height, width ) + width*.1*(rand-.5) +height*.1i*(rand-.5);
        end
        
        height=max(imag(rect))-min(imag(rect));
        width=max(real(rect))-min(real(rect));
        if info
            %fprintf('\n%e',max(height,width)); 
            %commented out, this is boring
        end
        
    end
    
    %at this point, we will have either one* zero or one* singularity.
    
    %polish things up with this function:
    [zeros, sings, zeroCount, singCount] = ZeroSingPolish( mean(rect), f, df, thresh, N^2 );
    
%     countOut=countZerosSingsRect( f, df, rect, N, 2 );
%     count=countOut(1);
%     count_X_z_0=countOut(2);
%     if count>0 %it's a zero
%         zeros=count_X_z_0/count;
%         zeroCount=count;
%     else %it's a singularity
%         sings=count_X_z_0/count;
%         singCount=count;
%     end
    
end