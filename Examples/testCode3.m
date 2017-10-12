%to do:
%- make oscillatory quadrature routine
%- remove unused matlab files
%- neaten up the kernel calls, decide on a framework for this
%- broaden incident waves, especially to beam-source type

clear classes;
addpath CompGaussFiles;
%wavenumber
kwave=5;

%define the screen
vertices=[0 0       %first vertex
          2*pi 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%inident plane wave
uinc=planeWave(kwave,[1 1]./sqrt(2));
%uinc=pointSource(kwave,[pi/2 -1]);

%with RHS data
f=DirichletData(uinc,Gamma);

%make an HNA basis on Gamma
pMax=4; nLayers=2*pMax; sigmaGrad=0.15;
V=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);

%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

%create integrator ojcect to do Galerkin integrals:
CG=CompGauss(kwave,20);

%initialise Galerkin matrix:
GalerkinMatrix=zeros(V.numEls);
GalerkinRHS=zeros(V.numEls,1);

basisSize=V.numEls; %copy this to avoid communication overhead
for n=1:basisSize
    
    for m=1:V.numEls
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
        GalerkinMatrix(n,m)=CG.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
    GalerkinRHSIntegrals2{n}=L2(S*GOA,V.el(n));
    GalerkinRHS(n)=CG.eval(GalerkinRHSIntegrals1{n})-CG.eval(GalerkinRHSIntegrals2{n});
    fprintf('%d of %d rows complete\n',n,basisSize);
end

coeffs=GalerkinMatrix\GalerkinRHS;
v_N=Projection(coeffs,V);
s=linspace(0,Gamma.L,100*kwave);
%plot(s,v_N.eval(s)+GOA.eval(s)); ylim([-3*kwave 3*kwave]); xlim(Gamma.supp);
plot(s,v_N.eval(s)+GOA.eval(s)); xlim(Gamma.supp); ylim([-3*kwave 3*kwave]); 