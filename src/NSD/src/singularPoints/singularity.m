classdef singularity
    %information about singularity
    
    properties
        position   %can be anonymous function, or just a vector,
                    %depending on the dimension of the singularity
        blowUpType    %can be 'log', or 'a' where x^a is the blowup type
         %nearlySingular  = false %true if function behaves nearly singularly
%         singDist_a = [] %distances of singularity from endpoint of domain
%         singDist_b = [] %distances of singularity from endpoint of domain
        distFun %  anonymous function which returns distance from singularity, given t
    end
    
    methods
        function self = singularity(position, blowUpType, distFun)
%             if nargin > 2
%                 self.nearlySingular = true;
%                 self.singDist_a = singDist_a;
%                 self.singDist_b = singDist_b;
%             end
            self.position=position;
            self.blowUpType=blowUpType;
            if nargin == 2
                self.distFun = @(r) abs(r-position);
            else
                self.distFun=distFun;
            end
        end
    end
    
end

