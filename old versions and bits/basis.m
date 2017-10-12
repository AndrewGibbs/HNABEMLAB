classdef (Abstract) basis <handle
    %have children HNA and hp
    
    properties
        numEls
        el
        side
    end
    
    methods
        function V=plus(V1,V2)
            %should introduce a check here for matching elements, to be
            %remvoed
            V.numEls=V1.numEls+V2.numEls;
            V.els=[V1.els V2.els];
            if isequal(V1.side,V2.side)
                V.side=V1.side;
            else
                %with this approach, understanding which basis element
                %points to which side is lost.
                V.side={V1.side,V2.side};
            end
        end
        
        %can't overload indexing function, as this already has a meaning
%         function out=subsindex(self,n)
%             out=self.el(n);
%         end
    end
    
end

