%have changed single mesh and hpStandard mesh files, to include pMax on
%each mesh element. Check this doesn't mess stuff up.

clear classes;
addpath ..;
addPaths();
%wavenumber
kwave=1000;

%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%inident plane wave
uinc=planeWave(kwave,[1 1]./sqrt(2));

%with RHS data
f=DirichletData(uinc,Gamma);

%make an HNA basis on Gamma
pMax=4; nLayers=2*(pMax+1); sigmaGrad=0.15; alphaDist=2;
VHNA=HNAsingleMesh(Gamma,pMax,kwave,alphaDist, nLayers, sigmaGrad);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

X = getColPoints( VHNA );
%X = ChebyshevRoots( DOFs, 'Tn', [Gamma.supp(1) Gamma.supp(2)] ).';

preColMatrix=[];
elCount=1;

%create collocation matrix
for b=VHNA.el
    Sb=S*b;
    preColMatrix=[preColMatrix Sb.col(X)];
    elCount=elCount+1;
end

%create RHSa, no integrals here
%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);% so can be computed exactly
RHSa=f.eval(X);

%create RHSb, these integrals will also have to be solved
Sgoa=S*GOA;
preColRHSb=Sgoa.col(X);

% initialise integration solvers
%CG=CompGaussBasic(kwave,1000,1000);
NSD=NSDlinearPhase(10,15);

%ColMatrixCG=zeros(size(preColMatrix));
ColMatrixNSD=zeros(size(preColMatrix));
ColRHSb=zeros(DOFs,1); ColRHSbGC=ColRHSb; ColRHSbNSD=ColRHSb;
for n=1:DOFs
%    CG_=CG;
    NSD_=NSD;
    fprintf('Computing column %d of %d\n',n,DOFs);
    for m=1:DOFs
        %ColMatrixCG(m,n)=CG_.eval(preColMatrix(m,n));
        ColMatrixNSD(m,n)=NSD_.eval(preColMatrix(m,n));
    end
    %ColRHSbGC(n)=CG_.eval(preColRHSb(n));
    ColRHSbNSD(n)=NSD_.eval(preColRHSb(n));
end
% 
% RelErr=abs(ColMatrixCG-ColMatrixNSD)./abs(ColMatrixCG);
% err=abs(ColMatrixCG-ColMatrixNSD);
% plot(abs(ColRHSbGC-ColRHSbNSD),'x');
%imagesc(err);
%mix the RHSs, f(x)-A\Psi(x):
ColRHS=RHSa-ColRHSbNSD;
coeffs=ColMatrixNSD\ColRHS;
v_N=Projection(coeffs,VHNA);
s=linspace(0,Gamma.L,100*kwave);
%plot(s,v_N.eval(s)+GOA.eval(s)); ylim([-3*kwave 3*kwave]); xlim(Gamma.supp);
plot(s,v_N.eval(s)+GOA.eval(s)); xlim(Gamma.supp); ylim([-3*kwave 3*kwave]); 