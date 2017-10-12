classdef BEMintegral1D < BEMintegral %again, not sure if this should be abstract or not
    
    properties
        ColPoint
        SubIntegrals
        type='mixed'
    end
    
    methods
        function self=BEMintegral1D(kernel, funciton, point) %Af(x) type situation
            self.ColPoint=point;
            self.domain=baseFn.supp;
            
            %this line will cause errors:
            BEMintegral1D.oscillator=self.basFn.oscillator*kernel.oscillator;
            
            if point.side == baseFn.side
                %should also have something in here for near singularities
                
                %if ColPoint is inside integration domain
                if point>self.domain(1) && point<self.domain(2)
                    self.type='SingularInsideDomain';
                    
                    %now split integral into two components, singular at
                    %edges
                    basFnL=copy(funciton);
                    basFnL.domain=[funciton.supp(1) point];
                    self.SubIntegrals{1}=BEMintegral1D(basFnL, point);
                    basFnR.domain=[point funciton.supp(2)];
                    self.SubIntegrals{2}=BEMintegral1D(basFnR, point);
                elseif point==self.domain(1)
                    self.type='SingularLeftDomain';
                elseif point==self.domain(2)
                    self.type='SingularRightDomain';
                else
                    error('Cannot classify integral');
                end
                
                %at this point, split up into subintegrals which can be
                %done nicely:
            else
                self.type='Smooth';
            end
        end
        
        function val=eval()
            %need to write some integration bits
        end
           
    end
    
end

