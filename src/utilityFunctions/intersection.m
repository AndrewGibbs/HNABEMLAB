function out = intersection( A,B )
    %ORIGINAL intersect code, returns 1 if AnB non-empty, zero otherwise.
    %boils down to four cases
%     yn=0;
% 
%     if A(1)<=B(1) && B(1)<A(2)
%         yn=1;
%     end
%     
%     if B(1)<=A(1) && A(1)<B(2)
%         yn=1;
%     end
%     
%     if B(1)<=A(1) && A(2)<=B(2)
%         yn=1;
%     end
%     
%     if A(1)<=B(1) && B(1)<A(2)
%         yn=1;
%     end
% 
%     if A(1)>=A(2) || B(1)>=B(2)
%         yn=0;
%    end
    z=min(A(2),B(2));
    y=max(A(1),B(1));
    
    if y>z
        out=[];
    elseif y==z
        out=y;
    else
        out=[y z];
    end

end