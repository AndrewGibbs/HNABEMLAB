%changes to make:
%-%decide on a convention for domain ordering, and adjust code accordingly
%-%should put a check in SingleLayerR2 to make sure function is on right
%domain, similar with L2 inner products

clear classes;
addpath CompGaussFiles;
kwave=10;

%define the first screen
vertices=[0 0       %first vertex
          2*pi 0];     %second vertex
%define the second screen
vertices2=[0 -1       %first vertex
          2*pi -1];     %second vertex
      
%create 'edge' object for the screen
Screen1=edge(vertices);
Screen2=edge(vertices2); %shift it a bit for the second obstacle

%define a density on the second obstacle
g=density(Screen2,@(s) 1);

%beam incidence is defined as
beamInc=SingleLayerR2(kwave,g,Screen2);

%with RHS Dirichlet data
f=DirichletData(beamInc,Screen1);

%take GOA approx of beam source on Screen1
beamGOA=GeometricalOpticsApprox(beamInc,Screen1);

%make an standard hp basis on Gamma
pMax=4; nLayers=2*pMax; sigmaGrad=0.15; hMax=2*pi/(2*kwave);
Vhp=hpStandardBasis(Screen1,pMax, hMax, nLayers, sigmaGrad);

%define the single layer 'operator' object
S=SingleLayer(kwave,Screen1);

%construct integrator objects
%composite/generalised Gauss, for 2D integrals
CG2d=CompGauss(kwave,30);
beamGOA.integratorDef=CG2d;
%brute force composite Gauss, for 3D integrals
%CG3d=CompGaussBasic(kwave,1000,30);

%initialise Galerkin matrix:
GalerkinMatrix=zeros(Vhp.numEls);
GalerkinRHS=zeros(Vhp.numEls,1);
basEls=Vhp.numEls;

for n=1:basEls
    for m=1:basEls
        GalerkinMatrixIntegrals{n,m}=L2(S*Vhp.el(m),Vhp.el(n));
    end
    GalerkinRHSIntegrals1{n}=L2(f,Vhp.el(n));
    %GalerkinRHSIntegrals2{n}=L2(S*beamGOA,V.el(n));
end

for n=1:basEls
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
coeffs=GalerkinMatrix\GalerkinRHS;
v_N=Vhp.project(coeffs);
% figure(1);
% plot(s,v_N.eval(s)); xlim(Screen1.supp); ylim([-.2*kwave .2*kwave]); 

%now construct an HNA basis
VHNA=HNAoverlappingMesh(Screen1,pMax, kwave, nLayers, sigmaGrad);
%create set of nodes with equal points on each mesh element
[s1,~]=CG2d.meshQuad(VHNA.mesh{1}); [s2,~]=CG2d.meshQuad(VHNA.mesh{2});
s=sort([s1; s2;]);
vHNA=VHNA.leastSquares(s,v_N.eval(s)-CG2d.eval(beamGOA,s));
plot(s,v_N.eval(s),'-',s,vHNA.eval(s)+CG2d.eval(beamGOA,s),':',s,CG2d.eval(beamGOA,s),'-.');
legend('hp approx','HNA least squares','POA');
ylim([-1 -0.3]); xlim([0 Screen1.L]);
%error plot of HNA approximating standard hp solution
%semilogy(s,abs(v_N.eval(s)-CG2d.eval(beamGOA,s)-vHNA.eval(s)))
