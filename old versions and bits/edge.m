classdef edge < handle
    %Screen object
    
    properties
        P1 %endpoint
        P2 %endpoint
        L   %length of screen
        dSv  %gradient of screen
        nv   %outward normal
    end
    
    methods
        
        %constructor
        function self = edge(vertices)
            if ~isequal(size(vertices),[2,2])
                error('Input for screen must be two by two sets of vertices');
            end
            %store endpoints
            self.P1=vertices(1,:);
            self.P2=vertices(2,:);
            self.L=norm(self.P1-self.P2);
            self.dSv=(self.P2-self.P1)/self.L;
            self.nv=[self.dSv(2) -self.dSv(1)];
        end
        
        function y=trace(self,s)
            s=s(:); %convert to upright vector
            if min(size(s))~=1
                error('Trace input must be a vector');
            end
            if ~(min(s)<0 || max(s)>self.L)
                y=repmat(self.P1,size(s))+s*self.dSv;
            else
                error('Parameter value outside of range');
            end
        end
        
        function y=normal(self,s)
            s=s(:); %convert to upright vector
            if min(size(s))~=1
                error('Normal trace input must be a vector');
            end
            if ~(min(s)<0 || max(s)>self.L)
                y=repmat(self.nv,size(s));
            else
                error('Parameter value outside of range');
            end
        end
        
    end
    
end