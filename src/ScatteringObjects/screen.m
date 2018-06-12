classdef screen < polygon
    %screen, which is like a one-sided polygon
    
    properties
        %side
        cumL
    end
    
    methods
        function self = screen(vertices)
           %self.side{1} = edge(vertices);
           self = polygon(vertices);
           self.cumL = self.side{1}.L;
        end
    end
end

