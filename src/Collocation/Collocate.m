classdef Collocate < handle
    %class function to simplify the collocation process
    
    properties
        sampleRatio
        ruleType
        pt
        edgeCol
    end
    
    methods
        function self = Collocate(Vbasis, sampleRatio, rule)
    
    
            self.pt = collocationPoint(blankMesh(),0,[],1,1,false);

             if ~isa(Vbasis.obstacle,'edge')
                 MultiCollocation(@(x) Collocate(x, sampleRatio, rule), Vbasis, self);
                 return;
             end
             
             %now assume collocating on a single side
             
             if sampleRatio < 1
                 error('Oversampling param must be at least one');
             end
             self.sampleRatio = sampleRatio;
             
             %decide on which quad rule to use
             if isequal(rule,'C')
                 self.ruleType = 'Chebyshev';
                 rule = @(pts) self.ChebRule(pts);
             elseif iseequal(rule,'U')
                 self.ruleType = 'Uniform';
                 rule = @(pts) self.UniformRule(pts)
             else
                 self.ruleType = 'User defined';
             end
             
             if isa(Vbasis,'HNAoverlappingMesh')
                  [E_, M] = Vbasis.mimicSingleMesh;
                  E = E_.el;
                  overlapElFlag = true;
             else
                E=Vbasis.mesh.el;
                M=Vbasis.meshDOFs;
                overlapElFlag = false;
             end
             
            sCount = 0;
            for m=1:length(M)
                pts=ceil(M(m)*sampleRatio);
                [s,w] = rule(pts);
             
                sSubCount = 0;
                for s_=s.'
                    sCount = sCount + 1;
                    sSubCount = sSubCount+1;
                    self.pt(sCount) = collocationPoint( E(m), .5*(s_+1), [], m, E(m).width*w(sSubCount)*.5, overlapElFlag);
                end
             end
        end
        
    end
    
    methods (Static)
        function [x,w] = ChebRule(pts)
            x=sort(cos(pi*(2*(1:pts)-1)/(2*pts))).';
            w = RiemannWeights(x,-1,1);
        end
        
        function [s,w] = UniformRule(pts)
            s=linspace(-1,1,pts+2).'; %add two extra points (endpoints)
            s=s(2:(end-1)); %delete the extra two points
            w = ones(pts,1)/pts;
        end
    end
end

