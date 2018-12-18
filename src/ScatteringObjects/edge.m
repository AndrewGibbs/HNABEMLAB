classdef edge < scatteringObject
    %Screen object
    
    properties
        P1 %endpoint
        P2 %endpoint
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
            self.supp=[0 self.L];
            self.numSides = 1; % screen has one side by default
        end
        
        function y=trace(self,s,~)
            if min(size(s))~=1
                error('Trace input must be a vector');
            end
            [sizeL sizeR] = size(s);
            s=s(:); %convert to upright vector
%             if ~(min(s)<0 || max(s)>self.L)
                 y=repmat(self.P1,size(s))+s*self.dSv;
%             else
%                error('Parameter value outside of range');
%             end
            if sizeL < sizeR %vector wasn't originally upright
                y=y.';
            end
        end
        
        function y=normal(self,s)
            s=s(:); %convert to upright vector
            if min(size(s))~=1
                error('Normal trace input must be a vector');
            end
%             if ~(min(s)<0 || max(s)>self.L)
                 y=repmat(self.nv,size(s));
%             else
%               error('Parameter value outside of range');
%             end
        end
        
        function R = dist(self,s,t)
            %returns distance between two points
            %** should eventually be moved into a superclass 'obstacle',
%             s=s(:);
%             t=t(:);
            if ~isequal(size(s),size(t)) && max(size(t))>1 && max(size(t))>1
                error('param vectors for distance function must be compatible');
            end
            R = abs(self.trace(s)-self.trace(t));
        end
        
        function R = distAnal(self,s,t,deriv,sGEt,~,~)
            if ~isequal(size(s),size(t)) && (max(size(t))>1 && max(size(s))>1 )
                error('param vectors for distance function must be compatible');
            end
%             s=s(:);
%             t=t(:);
            % sGEt is true if s>=t
            if nargin == 4
%                 sgn_ = sign(real(s-t));
%                 if min(sgn_) ~= max(sgn_) || ismember(0,sgn_)
%                     error('cannot guess analytic extension of distance function at zero')
%                 else
%                     sgn = min(sgn_);
%                 end
                error('need to know if s><t, to choose correct analytic extension');
%             elseif sGEt
%                 sgn = - 1;
%             else
%                 sgn= 1;
             end
            
%these are partial derivatives w.r.t. t
            switch deriv
                case 0
                    if sGEt
                        R=s-t;
                    else
                        R=t-s;
                    end
                case 1
                    if sGEt
                        R=-1;
                    else
                        R=1;
                    end
                otherwise
                    R=0;
            end
            %R = {@(s,t) sgn*(s-t), @(s,t) sgn, @(s,t) 0};
        end
          
        function draw(self)
           s = linspace(0,self.L);
           Y = self.trace(s);
           Y1 = Y(1,:);
           Y2 = Y(2,:);
           plot(Y1,Y2,'k','LineWidth',3);
        end
        
%         function sourceCheck(self,source)
%             %code which ensures wave is in upper-half plane
%             if isa(source,'SingleLayerR2') %beam source basically
%                 dot1=(self.P1-source.domain{2}.P1)*(self.nv.');
%                 dot2=(self.P1-source.domain{2}.P2)*(self.nv.');
%                 if sign(dot1)==-1*sign(dot2)
%                     error('Source appears to cross half-line, havent coded for this');
%                 elseif sign(dot1)>0
%                     self.nv=-self.nv;
%                     warning('Normal direction changed to face source');
%                 end
%             elseif isa(source,'planeWave')
%                 if source.d*(self.nv.')>0
%                     self.nv=-self.nv;
%                 end
%             elseif isa(source,'pointSource')
%                 if (self.P1-source.s)*(self.nv.')>0
%                     self.nv=-self.nv;
%                 end
%             else
%                 error('Source not classified');
%             end
%         end
    end
    
end