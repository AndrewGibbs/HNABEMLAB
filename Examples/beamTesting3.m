%% Description of HNA BEM LAB v1:
%This code constructs two parallel screens, one with a given surface
%density, and computes the surface density on the other.

% This code has been written as generally as possible, in an object
% oriented (OO) framework.
% EXAMPLE 1 OF WHY THIS FRAMEWORK IS NICE:
% Any HNA code requires an incident field (uinc). In particular, the
% Dirichlet/Neumann traces, and GOA are required. Previous HNA code has
% been designed around a plane wave incidence, but this may be considered a
% CHILD CLASS of a more general class of incident fields. All other
% children, point source, beam source, Herglotz-type incidences, share
% similar properties, although these are defined differently. At the top
% level, this code just requires (for e.g.) uinc.GOA, hence the same code
% can solve for any incident field.

% EXAMPLE 2 OF WHY THIS FRAMEWORK IS NICE:
%A second example would be the 'basis' class,
%  which has children 'HNAsingleMesh', 'HNAoverlappingMesh' and
% 'hpStandardBasis'. The code is written for a general basis, and can loop
% over each element, so may be applied to each basis.

%As integration is a key component of HNA, there is an 'integrator' class,
%alongside a range of different integeral classes, designed to capture all
%of the information about an integral (kernel(s), domain, support,
%functions) before computing it, by feeding it to an integrator class.
%This has made debugging it really easy, and is actually quite neat.
%Currently there are two integrator classes, 'CompGauss' and
%'CompGaussBasic'. The latter is a brute force Composite Gaussian rule,
%which has a fixed number of points per wavelength. The former is similar,
%but uses more sophisticated Generalised Gauss + Duffy routines for singularities.

%there are loads of classes, regrettably some (or at least some of their 
%properties) are redundant now and don't %get used, and the framework is not yet perfect.
%One disadvantage to this general OO framework,
% if I change a general class definiiton to fix one
%piece of code, this typically creates an error/bug elsewhere. I would
%advise to just create new classes, which may be children of other classes.
%For example, a 'Filon' class, which (like 'CompGauss' and
%'CompGaussBasic') are children of the 'integrator' parent class. This
%should also account for oscillations in the integrals.

%two key classes are BoundaryFunction, and BoundaryIntegral, which do
%pretty much what they say on the tin. As an example, for point source
%incidence, the POA 
%is a BoundaryFunction, and for beam source, the POA is BoundaryIntegral,
%but the integrator class must be written to distinguish between these.
% BoundaryFunction has the method 'eval', which evaluates at a point s. To
% evaluate a BoundaryIntegral, you need to also pass an 'integrator' class,
% to do the integral (evaluation is a more ambiguous term in this case).
% The density for the beam source is a child of BoundaryFunction, and the
% projection onto the HNA space (possible either via galerkin or least
% squares) is also a child of BoundaryFunction. This should lend it's self
% well to iterative solutions at some point, and at the second iteration
% one must write something like
%g=density(Screen2,@(s) v_HNA1.eval(s)+beamGOA1.eval(s));


%% General changes to make:
%should put a check in SingleLayerR2 to make sure function is on right
%domain, similar with L2 inner products
%oscillatory quadrature
%need to account for normal not pointing towards source


%% setup
clear classes; %clears more stuff than just 'clear' command
addpath CompGaussFiles; %files for composite Gauss, from PhD
%user defined parameters:
%-----------------------------------------------------------------------%
kwave=5; %wavenumber
pMax=1; nLayers=2*pMax; sigmaGrad=0.15; %basis parameters
%define the screens
vertices1=[0 0       %first vertex
          2*pi 0];     %second vertex
      
vertices2=[0 -1       %first vertex
          2*pi -1];     %second vertex
 %have to define source to be below screen we solve for, as formulation
 %isn't quite correct yet, and normal direction is defined pointing
 %DOWNWARDS
      
%-----------------------------------------------------------------------%

%create 'edge' object for the screen
Screen1=edge(vertices1); %\Gamma in our notation
Screen2=edge(vertices2); %\gamma in our notation
%(this edge class was intended to at somepoint be generalised to polygons)

%define an instance of density class on the second obstacle
g=density(Screen2,@(s) 1); %\varphi in our notation

%beam incidence is defined as
beamInc=SingleLayerR2(kwave,g,Screen2);
%SingleLayerR2 is a child of waveR2 class, which basically is a
%generalisation of the incident fields. At this stage, the incident field
%beamInc lives is defined on R2, this can be observed by checking 
%beamInc.domain, which contains the information that it maps from R2 to an
%edge.

%with RHS Dirichlet data
f=DirichletData(beamInc,Screen1);
%Dirichlet trace data only needed for the screen. A similar check to
%before, typing f.domain{1} and f.domain{2}, shows where f (which is an integral)
%maps to and from respectively. For nested integrals, domain{1} contains
%the domain onto which the BIO AB...Zf(x) maps, i.e. where x lives.

%take GOA approx of beam source on Screen1
beamGOA=GeometricalOpticsApprox(beamInc,Screen1);
%essentially twice the Neumann data here

%make an HNA basis on screen1
V=HNAoverlappingMesh(Screen1,pMax,kwave, nLayers, sigmaGrad);
%V is an instance of HNAoverlappingMesh class, which is a child of the
%basis class.

%define the single layer 'operator' object
S=SingleLayer(kwave,Screen1);
%these operator objects may seem pretty abstract, but simple, as there's not
%much information stored inside of them. Just type 'S' into the command
%line to see this. The complicated bit, which I'm pretty proud of, is that
%you can multiple an operator S by a function f, S*f is a defined
%operation. Also if you have a BoundaryIntegral Af, then S*Af is defined
%too. I;m not sure how robust this is...

%% Galerkin projection
%construct integrator object
%composite/generalised Gauss, for 2D integrals
CG=CompGauss(kwave,12);
%the second argument denotes number of points per wavelength.

%initialise Galerkin matrix with a load of zeros:
GalerkinMatrix=zeros(V.numEls);
GalerkinRHS=zeros(V.numEls,1);
basEls=V.numEls; %number of basis elements of V

for n=1:basEls
    for m=1:basEls
        %create array of integrals, but do not evaluate them yet. This may
        %be interpreted as semi-symbolic, where the integrals all remain
        %symbolic, until the tools are provided to evaluate them.
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
        %V.el(n) corresponds to the n'th basis element of the basis V.
    end
    %create two RHS arrays of integrals:
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
    GalerkinRHSIntegrals2{n}=L2(S*beamGOA,V.el(n));
    %notice how S*beamGOA is defined, i.e. an operator can be applied to
    %another operator. This is pretty cool. What is also quite neat, is
    %that the L2 inner product is defined, for 'L2(f,g)', or 'L2(Af,g)'.
end

%now we compute the integrals that we have just constructed, using CG, our
%compositeGauss.
%I have set this up so the following loop can be run as parfor,
% this seems to either speeds things up quite a bit, or crashes it completely.
parfor n=1:basEls %loop over all basis elements
    CG_=CG;    %need to make a copy of integrator for parfor to work:
    for m=1:basEls
        GalerkinMatrix(n,m)=CG_.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHS(n)=CG_.eval(GalerkinRHSIntegrals1{n})-CG_.eval(GalerkinRHSIntegrals2{n});
    fprintf('row %d (of %d total rows) complete\n',n,basEls);
end
%solve system to get coefficients
coeffsHNA=GalerkinMatrix\GalerkinRHS;
%create a BoundaryFunction by projecting onto basis V with these
%coeffiecient (this can be done with any basis, and least squares is coded
%too):
v_HNA=V.project(coeffsHNA);
%now we add an integrator object to beamGOA, such that it can be evaluated
%like a function:
beamGOA.integratorDef=CG; %use same CompGauss as before
%now get points to evaluate solution at. Can pass the mesh of the basis to
%the CG routine (but not all integrator classes) to generate nodes on each
%mesh element
[s1,~]=CG.meshQuad(V.mesh{1}); [s2,~]=CG.meshQuad(V.mesh{2});
%the above line is overkill for a simple plot, but would be necessary for a
%least-squares projection.
s=sort([s1; s2;]);
%plot the solution
plot(s,real(v_HNA.eval(s)+beamGOA.eval(s)),s,imag(v_HNA.eval(s)+beamGOA.eval(s)));
xlim(Screen1.supp); ylim([-1 1]); 
%----------------------------------------------------------
%now compare coeffs against a least squares approx of solution
%add extra param:
hMax=2*pi/(2*kwave);
%now make standard basis:L
Vhp=hpStandardBasis(Screen1,pMax, hMax, nLayers, sigmaGrad);
basEls=Vhp.numEls;
%reset
GalerkinMatrixHNA=GalerkinMatrix; %save this
GalerkinMatrix=zeros(V.numEls);
GalerkinRHS=zeros(V.numEls,1);
for n=1:basEls
    for m=1:basEls
        GalerkinMatrixIntegrals{n,m}=L2(S*Vhp.el(m),Vhp.el(n));
    end
    GalerkinRHSIntegrals1{n}=L2(f,Vhp.el(n));
    %GalerkinRHSIntegrals2{n}=L2(S*beamGOA,V.el(n));
end

parfor n=1:basEls
    CG2d_=CG2d;
    %CG3d_=CG3d;
    for m=1:basEls
        GalerkinMatrix(n,m)=CG2d_.eval(GalerkinMatrixIntegrals{n,m});
%        GalerkinMatrix_(n,m)=CG3d_.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHS(n)=CG2d_.eval(GalerkinRHSIntegrals1{n});%-CG3d_.eval(GalerkinRHSIntegrals2{n});
    %GalerkinRHS_(n)=CG3d_.eval(GalerkinRHSIntegrals1{n});
    fprintf('%d of %d rows complete\n',n,basEls);
end

coeffsStd=GalerkinMatrix\GalerkinRHS;
%hp discrete sol'n:
v_hp=Vhp.project(coeffsStd);

%create points to approximate least squares at:
[s1,~]=CG2d.meshQuad(V.mesh{1}); [s2,~]=CG2d.meshQuad(V.mesh{2});
s=sort([s1; s2;]);

[~,coeffsHNAls]=V.leastSquares(s,v_hp.eval(s)-CG2d.eval(beamGOA,s));
x=1:length(coeffsHNAls);
semilogy(x,abs(coeffsHNAls-coeffsHNA)./abs(coeffsHNAls),'x');
dataHNAls=GalerkinMatrixHNA*coeffsHNAls;
dataError=abs(GalerkinRHS-dataHNAls)./abs(dataHNAls);
semilogy(x,abs(GalerkinRHS-dataHNAls)./abs(dataHNAls),'x');
[~,maxErrorIndex]=max(dataError);