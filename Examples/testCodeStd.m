%bugs/fixes needed:
    %need to define a decent kernel

clear classes;
addpath CompGaussFiles;
%wavenumber
kwave=5;

%define the screen
vertices=[0 0       %first vertex
          2*pi 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%uinc=planeWave(kwave,[1 1]./sqrt(2));
uinc=pointSource(kwave,[pi/2 -1]);

%with RHS data
f=DirichletData(uinc,Gamma);

%make a standard basis on Gamma
pMax=4; nLayers=2*pMax; hMax=2*pi/(2*kwave); sigmaGrad=0.15;
V=hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);

%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

%create integrator ojcect to do Galerkin integrals:
CG=CompGaussBasic(kwave,200);
%CG=CompGauss(kwave,20);
%initialise Galerkin matrix:
GalerkinMatrix=zeros(V.numEls);
GalerkinRHS=zeros(V.numEls,1);

%create integrals to be solved
for n=1:V.numEls
    for m=1:V.numEls
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
    end
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
end

%now solve integrals in parallel
basisSize=V.numEls; %copy this to avoid communication overhead
for n=1:basisSize
    CGcopy=CG; %copy this to avoid communication overhead
    for m=1:basisSize
        GalerkinMatrix(n,m)=CGcopy.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHS(n)=CG.eval(GalerkinRHSIntegrals1{n});%-CG.eval(GalerkinRHSIntegrals2{m});
    fprintf('row %d complete \n',n);
end

%symmetry error test
% QuadErr=abs(GalerkinMatrix-(GalerkinMatrix.'))./abs(GalerkinMatrix);
% imagesc(QuadErr);
% [maxVal, maxEntry]=MatMax(QuadErr)

%construct solution:
coeffs=GalerkinMatrix\GalerkinRHS;
v_N=Projection(coeffs,V);
s=linspace(0,Gamma.L,1000);
plot(s,v_N.eval(s),s,GOA.eval(s)); ylim([-2*kwave 2*kwave]); xlim(Gamma.supp);
legend('v_N','POA');

% for n=1:V.numEls
%     for m=1:V.numEls
%         BasisProduct(n,m)=CG.eval(L2(V.el(m),V.el(n)));
%     end
% end
