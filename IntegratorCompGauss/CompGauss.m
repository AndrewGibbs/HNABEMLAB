classdef CompGauss < Integrator
    %Composite Gauss integrator, basic and inefficient, for testing
    
    properties
        wavelength   %will want to increase points per wavelength
        pointsPerWavelength=20
        GradingParam=0.15
        nearnesParam=0.2;   %max distance of two mesh elements to be classified as 'near'
        nLayers=12          %number of graded layers, for rare occasions when grading needed
        MemTol=1E5          %maximum number of points to sum at once
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
%                 if integral.toR2==0
%                    X=integral.domain{1}.trace(X); 
%                 end

                supp=integral.Fn.supp;

                %now create weights and nodes:
                [t, w]=gauss_quad_wave_split(supp(1), supp(2), self.wavelength, self.pointsPerWavelength );

%                 if integral.toR2==1
%                     Y_s=integral.domain.trace(t);
%                 else
                    Y_s=integral.domain{2}.trace(t);
                %end
                if integral.toR2==1

                    [X_1,Y_1]=meshgrid(X(:,1), Y_s(:,1));
                    [X_2,Y_2]=meshgrid(X(:,2), Y_s(:,2));

                    R=sqrt( (X_1-Y_1).^2 + (X_2-Y_2).^2 );
                    val=integral.evalR(R')*(w.*integral.Fn.eval(t));
                else
                    X=X(:);
                   val=(integral.kernel.evalFull(X,t,integral.domain{1},integral.domain{2},1)).'*(w.*integral.Fn.eval(t));
                end
            elseif isa(integral,'BEMintegral2D')
                if ~(isequal(integral.domain{1},integral.Fn1.domain) && isequal(integral.domain{2},integral.Fn2.domain))
                    error('Integral domains do not match function domains');
                end
                %this is where most of the action happens
                val=0;
                NearlySplitInts=self.splitIntoSubIntegrals(integral);
                %check for subintegrals of new set of integrals:
                Icount=1;
                for I_=NearlySplitInts
                    I=I_{1};%little bodge
                    if isequal(I.subIntegrals,[])
                        splitInts{Icount}=I;
                        Icount=Icount+1;
                    else
                        subSubInts=length(I.subIntegrals);
                        for j=1:subSubInts
                            splitInts{Icount}=I.subIntegrals{j};
                            Icount=Icount+1;
                        end
                    end
                end
                 X=[]; Y=[]; R=[]; W=[];
                 X_=[]; Y_=[]; R_=[]; W_=[];
%                 CGb=CompGaussBasic(5,500);
                %CG=CompGauss(5,20); 
                 Icount=1;
                 MemDumpCount=1;
                Icomp=zeros(size(splitInts));
                if isa(integral.Fn2,'PseudoBoundaryFunction')
                    hiddenWidth=integral.Fn2.BoundaryIntegral.domain{2}.L;
                    memAdj=self.pointsPerWavelength*ceil(hiddenWidth/self.wavelength);
                else
                    memAdj=1;
                end
                for I_=splitInts   %loop over all subintegrals
                    I=I_{1};%little bodge
                   switch I.type 
                       case 'diagSing'
                           [ x,y,r, w] = Duffy_quad_genG( I.supp{1}(1), I.supp{1}(2), 2*self.pointsPerWavelength, I.suppWidth(1));
                       case 'cornerSing'
                           if I.supp{1}(1)<I.supp{2}(1)
                              [x,y,r,w, bmx, ymb]=DuffyTouchLog(I.supp{1}(1),I.supp{1}(2),I.supp{2}(2),pi,self.pointsPerWavelength,self.GradingParam, self.nearnesParam, self.nLayers,I.suppWidth(1),I.suppWidth(2)); 
                           else
                              [y,x,r,w, ymb, bmx]=DuffyTouchLog(I.supp{2}(1),I.supp{2}(2),I.supp{1}(2),pi,self.pointsPerWavelength,self.GradingParam, self.nearnesParam, self.nLayers,I.suppWidth(2),I.suppWidth(1)); 
                           end
                       case 'nearSing'
                           b = geometricMiddle( I.supp{1}, I.supp{2} , I.suppWidth(1), I.suppWidth(2), self.nearnesParam);
                           if I.supp{1}(1)<I.supp{2}(1)
                              [x, y, r, w, bmx, ymb]= Graded_edges( I.supp{1}(1), I.supp{1}(2), b, I.supp{2}(1), I.supp{2}(2), pi, self.nLayers, self.GradingParam, self.pointsPerWavelength, I.suppWidth(1),I.suppWidth(2));
                           else
                              [y, x, r, w, ymb, bmx]= Graded_edges( I.supp{2}(1), I.supp{2}(2), b, I.supp{1}(1), I.supp{1}(2), pi, self.nLayers, self.GradingParam, self.pointsPerWavelength, I.suppWidth(2),I.suppWidth(1));
                           end                           
                       case 'smooth'
                           [ x,y,w]=Gauss_2D_quad_widths( I.supp{1}, I.supp{2}, self.pointsPerWavelength);
                           r=abs(x-y);
                       case 'unknown'
                           error('Integrals must be classified for CompGauss');
                   end
%                 if isequal(integral.domain{1},integral.domain{2}) %if both functions on same edge
%                     Icomp(Icount)=(w.')*(integral.kernel.eval(r).*I.Fn2.eval(y).*conj(I.Fn1.eval(x)));
%                 else
%                     Icomp(Icount)=(w.')*(integral.kernel.evalFull(x,y,integral.domain{1},integral.domain{2}).*I.Fn2.eval(y).*conj(I.Fn1.eval(x)));
%                 end
                   %Icomp(Icount)=w.'*(integral.kernel.eval(r).*I.Fn2.eval(y).*conj(I.Fn1.eval(x)));
                   %ICG(Icount)=CG.eval(I);
                    X=[X; x;]; Y=[Y; y;]; R=[R; r;]; W=[W; w;];
                    X_=[X_; x;]; Y_=[Y_; y;]; R_=[R_; r;]; W_=[W_; w;];
%                    Iest(Icount)=CGb.eval(I);
                    if length(X_)>self.MemTol/memAdj || Icount==length(splitInts)
                        %memory limit reached, or final integral reached -
                        %either way sum everything up
                        if isequal(integral.domain{1},integral.domain{2}) %if both functions on same edge
                            Icomp(MemDumpCount)=(W_.')*(integral.kernel.eval(R_).*integral.Fn2.eval(Y_).*conj(integral.Fn1.eval(X_)));
                        else
                            Icomp(MemDumpCount)=(W_.')*(integral.kernel.evalFull(X_,Y_,integral.domain{1},integral.domain{2}).*integral.Fn2.eval(Y_).*conj(integral.Fn1.eval(X_)));
                        end
                    MemDumpCount=MemDumpCount+1;
                    %reset weights and nodes
                     X_=[]; Y_=[]; R_=[]; W_=[];
                    end
                    Icount=Icount+1;
                end
                val=sum(Icomp);
%                 if isequal(integral.domain{1},integral.domain{2}) %if both functions on same edge
%                     val=(W.')*(integral.kernel.eval(R).*integral.Fn2.eval(Y).*conj(integral.Fn1.eval(X)));
%                 else
%                      val=(W.')*(integral.kernel.evalFull(X,Y,integral.domain{1},integral.domain{2}).*integral.Fn2.eval(Y).*conj(integral.Fn1.eval(X)));
%                 end
            elseif isa(integral,'InnerProduct1D')
                if length(integral.supp)==2 %check that integration domain has positive measure
                    [s, w]=gauss_quad_wave_split(integral.supp(1), integral.supp(2), self.wavelength, self.pointsPerWavelength );
                    val=w.'*(integral.Fn1.eval(s).*conj(integral.Fn2.eval(s)));
                else
                    val=0;
                end
            elseif isa(integral,'BEMintegral3D')
                if isequal(integral.domain{3},integral.domain{2})
                    error('This 3D integration bodge is designed for second and third domains to be disjoint');
                end
                %treat inner boundary integral as a function Gg(x)
                
                %3D integrals can be expensive. Might want to split it up.
                %estimate size of integral
%                 wavelengthsIn3Dintegral=ceil(integral.Fn1.suppWidth/self.wavelength)*ceil(integral.domain{2}.L/self.wavelength)*ceil(integral.Fn2.suppWidth/self.wavelength);
%                 MemEst=(wavelengthsIn3Dintegral*self.pointsPerWavelength)^3;
                 GgDom{1}=integral.domain{2}; GgDom{2}=integral.domain{3};
                 fDom{1}=integral.domain{1}; fDom{2}=integral.domain{2};
%                 %reconstruct BoundaryIntegral
                 Gg=BoundaryIntegral(integral.kernel{2}, integral.Fn2, GgDom);

%               %now fudge this so it acts like a function:
                Gg_as_Fn2=PseudoBoundaryFunction(Gg,self);
                integral2Dfudge=BEMintegral2D(fDom, integral.kernel{1}, integral.Fn1, Gg_as_Fn2);
                %send back into self, as 2D integral
                val=self.eval(integral2Dfudge);
            elseif isa(integral,'BEMintegral3D')
            
            elseif isa(integral,'BEMintegral1D')
                I=self.splitIntoSubIntegrals1D(integral);
                %......
            end

            
        end

        function Js=splitIntoSubIntegrals1D(self,integral)
            
                %check for subintegrals:
                if isequal(integral.subIntegrals,[])
                    Is=integral;
                    subInts=0;
                else
                    Is=integral.subIntegrals;
                    subInts=1;
                end
                
                
                nCount=1;
                for I=Is %may cause an error if can't loop over Is{n}
                    %determine how much integral needs to be split
                    if subInts==1
                        I_=I{1};
                    else
                        I_=I;
                    end
                    splits=ceil(I_.suppWidth/self.wavelength);
                    newWidth=I_.suppWidth/splits1;
                    endPoints1=linspace(I_.Fn.supp(1),I_.Fn.supp(2),splits1+1);
                    for n1=1:splits
                        %restrict functions to smaller subdomains
                        Gn=I_.Fn.restrictTo([endPoints1(n1),endPoints1(n1+1)]);
                        Js{nCount}=BEMintegral1D(integral.domain, I_.kernel, Gn);
                        nCount=nCount+1;
                    end
                end
        end
        
        function Js=splitIntoSubIntegrals(self,integral)
            
                %check for subintegrals:
                if isequal(integral.subIntegrals,[])
                    Is=integral;
                    subInts=0;
                else
                    Is=integral.subIntegrals;
                    subInts=1;
                end
                
                
                nCount=1;
                for I=Is %may cause an error if can't loop over Is{n}
                    %determine how much integral needs to be split
                    if subInts==1
                        I_=I{1};
                    else
                        I_=I;
                    end
                    splits1=ceil(I_.suppWidth(1)/self.wavelength);
                    newWidth1=I_.suppWidth(1)/splits1;
                    splits2=ceil(I_.suppWidth(2)/self.wavelength);
                    newWidth2=I_.suppWidth(2)/splits2;
                    endPoints1=linspace(I_.Fn1.supp(1),I_.Fn1.supp(2),splits1+1);
                    endPoints2=linspace(I_.Fn2.supp(1),I_.Fn2.supp(2),splits2+1);
                    for n1=1:(splits1)
                        for n2=1:(splits2)
                            %restrict functions to smaller subdomains
                            Gn1=I_.Fn1.restrictTo([endPoints1(n1),endPoints1(n1+1)]);
                            Gn2=I_.Fn2.restrictTo([endPoints2(n2),endPoints2(n2+1)]);
%                             Gn1=I_.Fn1.restrictTo([I_.Fn1.supp(1)+(n1-1)*newWidth1  I_.Fn1.supp(2)-(splits1-n1)*newWidth1]);
%                             Gn2=I_.Fn2.restrictTo([I_.Fn2.supp(1)+(n2-1)*newWidth2  I_.Fn2.supp(2)-(splits2-n2)*newWidth2]);
                            Js{nCount}=BEMintegral2D(integral.domain, I_.kernel, Gn1, Gn2);
                            nCount=nCount+1;
                        end
                    end
                end
        end

        
        function [s, w]=meshQuad(self,meshIn)
            %if a basis is passed by accident, just extract the mesh from
            %it:
            if isa(meshIn,'basis')
                meshIn=meshIn.mesh;
            end
            %account for double mesh cases:
            if iscell(meshIn)
                meshPoints=[];
                for m=1:length(meshIn)
                    meshPoints=[meshPoints meshIn{m}.points];
                    meshPoints=unique(sort(meshPoints));
                end
            else
                meshPoints=meshIn.points;
            end
                
            s=[]; w=[];
            for m=1:(length(meshPoints)-1)
                    [s_, w_]=gauss_quad_wave_split(meshPoints(m), meshPoints(m+1), self.wavelength, self.pointsPerWavelength );
                    s=[s; s_]; w=[w; w_];
%             for m=1:length(meshIn.el)
%                 [s_, w_]=gauss_quad_wave_split(meshIn.el(m).interval(1), meshIn.el(m).interval(2), self.wavelength, self.pointsPerWavelength, meshIn.el(m).width );
%                 s=[s; s_]; w=[w; w_];
%             end
            end
                
        end
    end
    
end