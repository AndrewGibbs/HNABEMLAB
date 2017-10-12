%have changed single mesh and hpStandard mesh files, to include pMax on
%each mesh element. Check this doesn't mess stuff up.

clear classes;
addpath ..;
addPaths();
%wavenumber
kwave=16;

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
pMax=5; nLayers=2*(pMax+1); sigmaGrad=0.15;
VHNA=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

X = getColPoints( VHNA, 0,1/sqrt(kwave));
COLs=length(X);

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

% initialise integration solver
NSD=NSDlinearPhase(20,20);

ColMatrixNSD=zeros(size(preColMatrix));
ColRHSb=zeros(COLs,1);
for m=1:COLs
    fprintf('Computing column %d of %d\n',m,COLs);
    for n=1:DOFs
        ColMatrixNSD(m,n)=NSD.eval(preColMatrix(m,n));
    end
    ColRHSb(m)=NSD.eval(preColRHSb(m));
end
fprintf('Condition number of matrix: %e',cond(ColMatrixNSD));
ColRHS=RHSa-ColRHSb;
coeffs=ColMatrixNSD\ColRHS;
v_N=Projection(coeffs,VHNA);
s=linspace(0,Gamma.L,100*kwave);
plot(s,v_N.eval(s)+GOA.eval(s)); xlim(Gamma.supp); ylim([-3*kwave 3*kwave]); 