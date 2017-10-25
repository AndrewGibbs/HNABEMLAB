classdef NSDv1 < Integrator
    %A numerical steepst descent integrator, which takes a 1D integral as
    %an input, and (provided the phase is linear), returns the value of the
    %integral, using numerical steepest descent.
    %assumes that the collocation point and the basis element are on the
    %same side, which is fine - if they aren't, then the the phase will not
    %be linear.
    
    properties
        Nsd %number of nodes to take along SD paths
        Ngauss  % number of nodes to take for standard Gaussian quadrature
        nLayers=12   %number of graded layers
        GradingParam=0.15 %geometric grading parameter for singularities
    end
    
    methods 
        function self=NSDlinearPhase(Nsd,Ngauss, nLayers,GradingParam)
            %set key parameters, if they are specified
            if nargin==0
                self.Nsd=10;
            else
                self.Nsd=Nsd;
            end
            
            if nargin>=2
                self.Ngauss=Ngauss;
            else
                self.Ngauss=self.Nsd;
            end
            
            if nargin>=3
                self.nLayers=nLayers;
            end
            if nargin>=4
                self.GradingParam=GradingParam;
            end
            
            %at this point might be worth generating one set of weights and
            %nodes for each type of integral required, and calling them 
        end
        
        function val=eval(self,I)
            %check input is an integral
            if ~isa(I,'Integral')
                error('NSD requires integrals to solve, input was not of Integral.m type');
            end
            %check there are no free variables left in integral
            if length(I.domain)~=I.dim
               error('size of domain not equal to number of dimensions of integral'); 
            end
            
            %test if integral is oscillatory:
            if isequal(I.oscillator,[])
                osc = zeros(I.dim,1);
            else
                %perhaps this could be a little more sophisticated:
                osc = I.width > 2*pi/I.oscillator.freq  ;
            end
            
            
            if I.dim==1
                if osc==1 %oscillatory integral, transform it into something that isn't
                    %can do 1D steepest descent easily
                    if I.analytic==true && strcmp(I.oscillator.phaseType,'linear')
                        J=self.NSDtransformLinearPhase(I);
                    end
                    
                    
                    
                else %split integral up s.t. singularities are at the corners
            end
            
            % J is a vector containing multiple integrals, each of which
            % are now non-oscillatory and exponentially decaying.
                
                %--------------integral is NOT oscillatory-------------
                width=I.suppWidth;
                if strcmp(I.type , 'Singular left')
                    [x,w,r]=genGaussLog( self.Ngauss, a, b, width, 'L' );
                elseif strcmp(I.type , 'Singular right')
                    [x,w,r]=genGaussLog( self.Ngauss, a, b, width, 'R' );    

                    %for nearly singulars,
                    %could subtract off singular bit, using negative
                    %weights, but would need smooth extension of Fn
                 elseif strcmp(I.type , 'Nearly singular left')
%                         [x1,w1,r1]=genGaussLog( self.Ngauss, I.ColPoint, b, b-I.ColPoint, 'L' );  
%                         [x2,w2,r2]=genGaussLog( self.Ngauss, I.ColPoint, a, a-I.ColPoint, 'L' );
%                         x=[x1; x2;]; w=[w1; -w2;]; r=[r1; r2;];
                      singR=a-I.ColPoint;
                      [r0, w0] = GradedQuad( self.Ngauss, self.nLayers, self.GradingParam );
                      r=singR+width*r0;	x=I.ColPoint+r;  w=w0*width;
                elseif strcmp(I.type , 'Nearly singular right')
%                         [x1,w1,r1]=genGaussLog( self.Ngauss, a, I.ColPoint, I.ColPoint-a, 'R' );
%                         [x2,w2,r2]=genGaussLog( self.Ngauss, b, I.ColPoint, I.ColPoint-b, 'R' );
%                         x=[x1; x2;]; w=[w1; -w2;]; r=[r1; r2;];  
                      singR=I.ColPoint-b; 
                      [r0, w0] = GradedQuad( self.Ngauss, self.nLayers, self.GradingParam );
                      r=singR + width*flipud(r0(:));   w=width*flipud(w0(:));
                      x=I.ColPoint-r;
                      %r=width*r0;	x=a+r;  w=w0*width;
                elseif strcmp(I.type , 'Smooth')
                    [x,w]=gauss_quad(a,b,self.Ngauss);
                    r=abs(x-I.ColPoint);
                else
                    error('Integral classification not recognised');
                end
                val=w.'*(I.kernel.eval(r).*I.Fn.eval(x));
            end
    end
    
    methods(Static)
                function NSDtransformLinearPhase(I)
                %can exploit linear-ness of phase.
                %make two copies of current integral
                a=I.domain{1}(1); b=I.domain{1}(2);
                
                    cCount=0;
                    newFreq=I.oscillator.freq*I.oscillator.phaseData(1);
                    for c=[a b]
                        %deform onto steepest descent paths
                        cCount=cCount+1;
                        Insd(cCount)=I;
                        Insd(cCount).oscillator=[];
                        Insd(cCount).domain=[0  inf];
                        if isequal(c,a)
                            pm=1;
                        else
                            pm=-1;
                        end
                        Insd(cCount).integrand = @(x) pm*exp(1i*newFreq*I.oscillator.phase(c)+I.oscillator.data(2)*I.oscillator.freq)/newFreq  *  I.integrandNonOsc(x/newFreq);
                        Insd(cCount).width=inf;
                        Insd(cCount).weightType='Gauss–Laguerre';
                    end 
                end
        
    end
    
end

