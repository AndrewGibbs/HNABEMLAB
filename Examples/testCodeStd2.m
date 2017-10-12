%bugs/fixes needed:
    %need to define a decent kernel

clear classes;
addpath Quadrature;
%wavenumber
kwave=5;

%define the screen
vertices=[0 0       %first vertex
          2*pi 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%uinc=planeWave(kwave,[1 1]./sqrt(2));
uinc=pointSource(kwave,[1/2 -1]);

%with RHS data
f=DirichletData(uinc,Gamma);

%make a standard basis on Gamma
pMax=2; nLayers=2*pMax; hMax=2*pi/(2*kwave); sigmaGrad=0.15;
V=hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);

%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

%create integrator ojcect to do Galerkin integrals:
%CG_=CompGaussBasic(kwave,1000);
CG=CompGauss(kwave,30);
%initialise Galerkin matrix:
GalerkinMatrixBad=zeros(V.numEls); GalerkinMatrixGood=GalerkinMatrixBad;
GalerkinRHSbad=zeros(V.numEls,1); GalerkinRHSgood=GalerkinRHSbad;
for n=1:V.numEls
    for m=1:V.numEls
        GalerkinMatrixIntegrals{n,m}=L2(S*V.el(m),V.el(n));
        GalerkinMatrixBad(n,m)=CG.eval(GalerkinMatrixIntegrals{n,m});
        %GalerkinMatrixGood(n,m)=CG_.eval(GalerkinMatrixIntegrals{n,m});
    end
    GalerkinRHSIntegrals1{n}=L2(f,V.el(n));
    GalerkinRHSbad(n)=CG.eval(GalerkinRHSIntegrals1{n});
    %GalerkinRHSgood(n)=CG_.eval(GalerkinRHSIntegrals1{n});
    fprintf('%d / %d\n',n,V.numEls);
end
% %symmetry error test
% QuadErr=abs(GalerkinMatrixBad-(GalerkinMatrixBad.'))./abs(GalerkinMatrixBad);
% imagesc(QuadErr);
% [maxVal, maxEntry]=MatMax(QuadErr)

%construct solution:
coeffs=GalerkinMatrixBad\GalerkinRHSbad;
v_N=Projection(coeffs,V);
s=linspace(0,Gamma.L,1000);
plot(s,v_N.eval(s),'.',s,GOA.eval(s));% ylim([-2*kwave 2*kwave]);
legend('v_N','POA');

% for n=1:V.numEls
%     for m=1:V.numEls
%         BasisProduct(n,m)=CG.eval(L2(V.el(m),V.el(n)));
%     end
% end
