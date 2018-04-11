classdef planeWave < waveR2
    %plane wave incident field
    
    properties
        d
    end
    
    methods
        
        function self=planeWave(kwave,direction)
           self.kwave=kwave;
           self.phaseLinear=direction;
           self.d=direction(:).';
           self.phasePD{1,1} = @(X) X*(self.d.');
           self.phasePD{2,1} = @(X) repmat(self.d(1),[size(X,1) 1]);
           self.phasePD{1,2} = @(X) repmat(self.d(2),[size(X,1) 1]);
           self.phasePD{2,2} = @(X) repmat(0,[size(X,1) 1]);
        end

%         function g=phase(self,X)
%            g=X*self.d;
%         end
        
        function [valOsc, val]=DirTrace(self,s,boundary)
            %standard def of plane wave
            valOsc=exp(1i*self.kwave*boundary.trace(s)*self.d.');
            val=ones(length(s),1);
        end
        
        function [val, valNonOsc]=NeuTrace(self,s,boundary)
            %derivative of plane wave, product with normal vector
            valNonOsc=1i*self.kwave*boundary.normal(s)*self.d.';
            val=valNonOsc.*self.DirTrace(s,boundary);
        end
%         
%         function [valOsc, val]=POA(self,s,boundary)
%             %actually, just need to be consistent with normal direction
%             val=2*self.NeuTrace(self,s,boundary);
%         end
    end
    
end