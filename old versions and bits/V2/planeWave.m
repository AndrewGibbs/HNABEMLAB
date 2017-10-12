classdef planeWave < waveR2
    %plane wave incident field
    
    properties
        d
    end
    
    methods
        
        function self=planeWave(kwave,direction)
           self.kwave=kwave;
           self.d=direction;
           self.oscillator=@(x) 1i*(x(1)*direction(1) + x(2)*direction(2));
        end
        
        function [valOsc, val]=DirTrace(self,s,boundary)
            %standard def of plane wave
            valOsc=exp(1i*self.kwave*boundary.trace(s)*self.d.');
            val=ones(length(s),1);
        end
        
        function [valOsc, val]=NeuTrace(self,s,boundary)
            %derivative of plane wave, product with normal vector
            val=1i*self.kwave*boundary.normal(s)*self.d.';
            valOsc=val.*self.DirTrace(s,boundary);
        end
%         
%         function [valOsc, val]=POA(self,s,boundary)
%             %actually, just need to be consistent with normal direction
%             val=2*self.NeuTrace(self,s,boundary);
%         end
    end
    
end