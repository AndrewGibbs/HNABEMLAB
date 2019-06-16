classdef collocationPoint
    % all the info you could ever want about a single point
    
    properties
        x    % in [0,L]
        meshIndex   % in |N
        distSideL
        distSideR
        distMeshL = []
        distMeshR = []
        side %should be mostly redudanct now...
        weight=1 % in (0,infty) option to include a weight for colloction point
        %overlapFlag
    end
    
    methods
        function self = collocationPoint(meshEl,s,side,meshIndex,weight,overlapFlag)
            %choose the collocation point, and all relevant relative
            %distances
            %self.overlapFlag = midElFlag;
            self.x = meshEl.interval(1)+meshEl.width*s;
            
            self.distMeshL = meshEl.width*s;
            self.distMeshR = meshEl.width*(1-s);
            self.distSideL = self.distMeshL + meshEl.distL;
            self.distSideR = self.distMeshR + meshEl.distR;
            if overlapFlag
                self.distMeshL = [];
                self.distMeshR = [];
            end
            self.side = side;
            self.meshIndex = meshIndex;
            %if nargin == 5
                self.weight = weight;
            %end
           
        end
    end
end

