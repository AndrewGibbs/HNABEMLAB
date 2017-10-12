%bugs:
%-- currently the CompGauss splits the integrals badly due to rounding
    %errors, and more integrals are computed than is needed.
%-- CompGauss doesn't give same result as CompGaussBasic for (S\phi,\Psi)
%-- Overlapping basis is constructed wrong, with elemnents oscillating in
    %both directions on both meshes

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
CGb=CompGaussBasic(kwave,500);
%initialise Galerkin matrix:
GalerkinMatrix=zeros(V.numEls); GalerkinMatrixb=GalerkinMatrix;
GalerkinRHS=zeros(V.numEls,1);  GalerkinRHSb=GalerkinRHS;

basisSize=V.numEls; %copy this to avoid communication overhead
for n=1:basisSize
    
    for m=1:V.numEls
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
    end
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
    GalerkinRHSIntegrals2{n}=L2(S*GOA,V.el(n));
end
for n=1:basisSize
    CGcopy=CG; %copy this to avoid communication overhead
    CGcopyb=CGb; %copy this to avoid communication overhead
    for m=1:basisSize
        %GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
        GalerkinMatrix(n,m)=CGcopy.eval(GalerkinMatrixIntegrals{n,m});
        GalerkinMatrixb(n,m)=CGcopyb.eval(GalerkinMatrixIntegrals{n,m});
    end
%     GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
%     GalerkinRHSIntegrals2{n}=L2(S*GOA,V.el(n));
    GalerkinRHS(n)=CGcopyb.eval(GalerkinRHSIntegrals1{n})-CGcopyb.eval(GalerkinRHSIntegrals2{n});
    %GalerkinRHSb(n)=CGcopyb.eval(GalerkinRHSIntegrals1{n})-CGcopyb.eval(GalerkinRHSIntegrals2{n});
    fprintf('row %d / %d complete \n',n,basisSize);
end

coeffs=GalerkinMatrixb\GalerkinRHS;
v_N=Projection(coeffs,V);
s=linspace(0,Gamma.L,100*kwave);
%plot(s,v_N.eval(s)+GOA.eval(s)); ylim([-3*kwave 3*kwave]); xlim(Gamma.supp);
plot(s,v_N.eval(s)); ylim([-3*kwave 3*kwave]); xlim(Gamma.supp);