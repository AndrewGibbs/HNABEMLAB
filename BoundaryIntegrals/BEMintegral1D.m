classdef BEMintegral1D < BEMintegral %again, not sure if this should be abstract or not
    
    properties
        ColPoint
        SubIntegrals
        phaseLinear       %[m c], phase function of oscillatory component exp(i*k*(mx+c))
        singDist=.2    %distance from singularity at which integral is considered singular
        supp
        suppWidth
        Fn
        osc
    end
    
    methods
        function self=BEMintegral1D(domain, kernel, Fn, point, integralType) %Af(x) type situation
            
            if nargin<=4
                self.type='unknown';
            else
                self.type=integralType;
            end
            self.kernel=kernel;
            self.ColPoint=point;
            self.domain=domain;
            self.supp=Fn.supp;
            self.suppWidth=Fn.suppWidth;
            self.Fn=Fn;
            wavelength=2*pi/self.kernel.kwave;
            
%             %now combine the phases:
%             phaseLinKer=kernel.colPhase(x,);
%             if isequal(phaseLinKer,[])
%                 self.phaseLinear=[];
%             else
%                 self.phaseLinear=phaseLinKer+[Fn.pm 0];
%             end
            
            
            %self.singDist=0.15*kernel.kwave;
            %singBall=[point - self.singDist  point + self.singDist];
            %self.phase=
            
            
            %if point.side == baseFn.side
                %should also have something in here for near singularities
                
                %if ColPoint is inside integration domain
                if point>self.supp(1) && point<self.supp(2)
                    self.type='Singular';
                    Fn1=Fn.restrictTo([self.supp(1) point]);
                    Fn2=Fn.restrictTo([point self.supp(2)]);
                    self.SubIntegrals{1}=BEMintegral1D(domain, kernel, Fn1, point, self.type);
                    self.SubIntegrals{2}=BEMintegral1D(domain, kernel, Fn2, point, self.type);
                elseif point == self.supp(1)
                    self.type='Singular left';
%                     Fn1=Fn.ResBasFn([self.supp(1) point]);
%                     self.SubIntegrals{1}=BEMintegral1D(domain, kernel, Fn1, point, integralType);
                elseif point == self.supp(2)
                    self.type='Singular right';
%                     Fn2=Fn.ResBasFn([point self.supp(2)]);
%                     self.SubIntegrals{2}=BEMintegral1D(domain, kernel, Fn2, point, integralType);
                elseif   Fn.supp(1)-self.singDist < point && point < Fn.supp(1)
                    self.type='Nearly singular left';
                elseif point < Fn.supp(2) + self.singDist && point > Fn.supp(2)
                    self.type='Nearly singular right';
                else
                    self.type='Smooth';
                end
                
                if ~strcmp(self.type,'Singular') && ~isequal(Fn.phaseLinear,[])
                    if point <= Fn.supp(1)
                        phaseLinKer=kernel.colPhaseLin(point,'xLy');
                    else
                        phaseLinKer=kernel.colPhaseLin(point,'yLx');
                    end
                    
                    %now combine the linear phases of kernel and Fn:
                    self.phaseLinear=phaseLinKer+Fn.phaseLinear;
                else
                    %phase is not a linear function, so leave it blank
                    self.phaseLinear=[];
                end
                
                
            %if integration domain is greater than one wavelength
            if isequal(self.phaseLinear,[])
                self.osc=[];
            else
                if self.suppWidth > (wavelength/abs(self.phaseLinear(1)))
                    self.osc=1;
                else
                    self.osc=0;
                end
            end
                
                
        end
           
    end
    
end

