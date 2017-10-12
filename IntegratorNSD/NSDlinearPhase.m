classdef NSDlinearPhase < Integrator
    %A numerical steepst descent integrator, which takes a 1D integral as
    %an input, and (provided the phase is linear), returns the value of the
    %integral, using numerical steepest descent.
    %assumes that the collocation point and the basis element are on the
    %same side, which is fine - if they aren't, then the the phase will not
    %be linear.
    
    properties
        Nsd=10 %number of nodes to take along SD paths
        Ngauss=15  % number of nodes to take for standard Gaussian quadrature
        
        %parameters for occasional graded quadrature:
        nLayers=12   %number of graded layers
        GradingParam=0.15 %geometric grading parameter
    end
    
    methods 
        function self=NSDlinearPhase(Nsd,Ngauss, nLayers,GradingParam)
            %set key parameters, if they are specified
            if nargin>=1
                self.Nsd=Nsd;
            end
            if nargin>=2
                self.Ngauss=Ngauss;
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
            if ~isequal(I.SubIntegrals,[])
                %integral has been split into subintegrals
                subInts=length(I.SubIntegrals);
                val=0;
                for j=1:subInts
                    val=val+self.eval(I.SubIntegrals{j});
                end
            else
                %no subintegrals, need to evaluate at this level
                if ~isa(I,'BEMintegral1D')
                    error('Can only handle 1D integrals');
                end
                if isequal(I.phaseLinear,[])
                    error('Phase must be linear, and included in integral object');
                end
                
                %endpoints of integral:
                a=I.supp(1); b=I.supp(2);
                if I.osc
                    %--------------integral is oscillatory-----------------
                    k=I.kernel.kwave;
                    m_g=I.phaseLinear(1);   c_g=I.phaseLinear(2);
                    if strcmp(I.type , 'Singular left')
                        [ x, w, r ] = NSDLogLinearPhase( self.Nsd, k,a,b,m_g,c_g,'L' );
                    elseif strcmp(I.type , 'Singular right')
                        [ x, w, r ] = NSDLogLinearPhase( self.Nsd, k,a,b,m_g,c_g,'R' );
                    elseif strcmp(I.type , 'Nearly singular left') || strcmp(I.type , 'Nearly singular right')
                        [ x, w, r ] = NSDnearLogLinearPhase( self.Ngauss, k, a, b, m_g, c_g, I.ColPoint, self.nLayers, self.GradingParam );
                    elseif strcmp(I.type , 'Smooth')
                        [ x, w ] = SDbasic( self.Nsd, k,a,b,m_g,c_g );
                    else
                        error('Integral classification not recognised');
                    end
                    %r being output is currently wrong, here's a quick fix:
                    if I.ColPoint<=I.supp(1)
                        r=x-I.ColPoint;
                    else
                        r=I.ColPoint-x;
                    end
                    [~,NonOscKer]=I.kernel.eval(r);
                    [~,FnNonOsc]=I.Fn.eval(x);
                    val=w.'*(NonOscKer.*FnNonOsc);
                else 
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
        end
    end
    
end

