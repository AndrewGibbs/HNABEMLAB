%changes to make:
%-%decide on a convention for domain ordering, and adjust code accordingly
%-%should put a check in SingleLayerR2 to make sure function is on right
%domain, similar with L2 inner products

clear classes;
addpath CompGaussFiles;
kwave=1;

%define the screen
vertices1=[0 0       %first vertex
          2*pi 0];     %second vertex
      
vertices2=[0 -1       %first vertex
          2*pi -1];     %second vertex
      
%create 'edge' object for the screen
Screen1=edge(vertices1);
Screen2=edge(vertices2); %shift it a bit for the second obstacle

%define a density on the second obstacle
g=density(Screen2,@(s) 1);

%beam incidence is defined as
beamInc=SingleLayerR2(kwave,g,Screen2);

%with RHS Dirichlet data
f=DirichletData(beamInc,Screen1);

%take GOA approx of beam source on Screen1
beamGOA=GeometricalOpticsApprox(beamInc,Screen1);

%make an HNA basis on Gamma
pMax=4; nLayers=2*pMax; sigmaGrad=0.15;
V=HNAoverlappingMesh(Screen1,pMax,kwave, nLayers, sigmaGrad);

%define the single layer 'operator' object
S=SingleLayer(kwave,Screen1);

%construct integrator objects
%composite/generalised Gauss, for 2D integrals
CG2d=CompGaussBasic(kwave,30);
%brute force composite Gauss, for 3D integrals
CG3d=CompGaussBasic(kwave,300,30);

%initialise Galerkin matrix:
GalerkinMatrix=zeros(V.numEls);
GalerkinRHS=zeros(V.numEls,1);
basEls=V.numEls;

for n=1:basEls
    for m=1:basEls
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
    end
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
    GalerkinRHSIntegrals2{n}=L2(S*beamGOA,V.el(n));
end

parfor n=1:basEls
    CG2d_=CG2d;
    CG3d_=CG3d;
    for m=1:basEls
        GalerkinMatrix(n,m)=CG2d_.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHS(n)=CG2d_.eval(GalerkinRHSIntegrals1{n})-CG3d_.eval(GalerkinRHSIntegrals2{n});
    fprintf('%d of %d rows complete\n',n,basEls);
end
coeffs=GalerkinMatrix\GalerkinRHS;
v_N=V.project(coeffs);
s=linspace(0,Screen1.L,100*kwave);
plot(s,v_N.eval(s)); xlim(Screen1.supp); ylim([-1 1]); 
