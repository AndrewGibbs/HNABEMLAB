classdef collocationPoint
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x    % in [0,L]
        meshIndex   % in |N
        distSideL
        distSideR
        distMeshL
        distMeshR
        side
    end
    
    methods
        function self = collocationPoint(meshEl,s,side,meshIndex)
            %choose the collocation point, and all relevant relative
            %distances
            
            self.x = meshEl.interval(1)+meshEl.width*s;
            self.distMeshL = meshEl.width*s;
            self.distMeshR = meshEl.width*(1-s);
            self.distSideL = self.distMeshL + meshEl.distL;
            self.distSideR = self.distMeshR + meshEl.distR;
            self.side = side;
            self.meshIndex = meshIndex;
           
        end
    end
end

