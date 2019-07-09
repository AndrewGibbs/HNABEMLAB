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
           self.phaseMaxStationaryPointOrder=0;
        end
        
        function uincR = getReflection(self,mirror)
            %get reflection of plane wave, as planeWave object
            errMess = '2nd input must be a PolygonalScatteringObject with 2 vertices';
            if ~isa(mirror,'Screen') && ~isa(mirror,'Edge') && ~isa(mirror,'MultiScreen')
                error(errMess);
            elseif ~isequal(size(mirror.vertices),[2 2])
                error(errMess);                
            end
            if mirror.nv*self.d.'<0
                n = mirror.nv;
            else
                n = -mirror.nv;
            end
            dR = self.d - 2*mirror.nv*(self.d*n.');
            uincR = planeWave(self.kwave,dR);
        end
        
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
        
        function g = phasePD(self,X,xDers,yDers)
            zeroOption = repmat(0,[size(X,1) 1]);
            switch xDers
                case 0
                    switch yDers
                        case 0
                            g = X*(self.d.');
                        case 1
                            g = repmat(self.d(2),[size(X,1) 1]);
                        otherwise
                            g = zeroOption;
                    end
                case 1
                        switch yDers
                            case 0
                                g = repmat(self.d(1),[size(X,1) 1]);
                            case 1
                                g = zeroOption;
                            otherwise
                                g = zeroOption;
                        end
                otherwise
                    g = zeroOption;
            end
        end
        
        function x = sourceVsnormal(self,side)
            x = sign((self.d)*(side.nv.'));
        end
        
        function val = eval(self,X1,X2)
            val = exp(1i*self.kwave*(self.d(1)*X1+self.d(2)*X2));
        end
    end
    
end