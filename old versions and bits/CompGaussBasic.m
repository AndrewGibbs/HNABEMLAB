classdef CompGauss < Integrator
    %Composite Gauss integrator, basic and inefficient, for testing
    
    properties
        wavelength   %will want to increase points per wavelength
        pointsPerWavelength
    end
    
    methods
        function self=CompGauss(kwave,pointsPerWavelength)
            if nargin==1
                self.pointsPerWavelength=10; %default
            else
                self.pointsPerWavelength=pointsPerWavelength;
            end
            self.wavelength=2*pi/kwave;            
        end
        
        function val=eval(self,integral, X)
            if isa(integral,'BoundaryIntegral')
                %ideally would look at properties of integral here...

                %check for if parametrised values have been given:
                if integral.toR2==0
                   X=integral.domain{end}.trace(X); 
                end

                supp=integral.boundaryFn.supp;

                %now create weights and nodes:
                [s, w]=gauss_quad_wave_split(supp(1), supp(2), self.wavelength, self.pointsPerWavelength );

                 Y_s=integral.domain.trace(s);

                [X_1,Y_1]=meshgrid(X(:,1), Y_s(:,1));
                [X_2,Y_2]=meshgrid(X(:,2), Y_s(:,2));

                R=sqrt( (X_1-Y_1).^2 + (X_2-Y_2).^2 );
                val=integral.kernel(R')*(w.*integral.boundaryFn.eval(s));
            elseif isa(integral,'BEMintegral2D') && isequal(integral.Fn1.domain,integral.Fn2.domain)
                supp_s=integral.Fn1.supp;
                supp_t=integral.Fn2.supp;
                [s, w_s]=gauss_quad_wave_split(supp_s(1), supp_s(2), self.wavelength, self.pointsPerWavelength, integral.Fn1.suppWidth );
                [t, w_t]=gauss_quad_wave_split(supp_t(1), supp_t(2), self.wavelength, self.pointsPerWavelength +1, integral.Fn2.suppWidth);
                [S,T] = kron2( s,t ); [wS,wT] = kron2( w_s,w_t );
                val=(wS.*wT).'*(integral.kernel(S,T).*integral.Fn1.eval(S).*conj(integral.Fn2.eval(T)));
            elseif isa(integral,'InnnerProduct1D')
                [s, w]=gauss_quad_wave_split(integral.supp(1), integral.supp(2), self.wavelength, self.pointsPerWavelength );
                val=w.'*(integral.Fn1.eval(s).*conj(integral.Fn2.eval(s)));
            end
        end
    end
    
end