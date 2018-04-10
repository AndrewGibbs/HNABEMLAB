clc;
clear classes;
addpath ..;
addPaths();
messageFlag=1;
%wavenumber
kwave=128;

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
pMax=4; nLayers=2*(pMax+1)-1; sigmaGrad=0.15; throwAwayParam=0;
VHNA=HNAsingleMesh(Gamma,pMax,kwave, throwAwayParam, nLayers, sigmaGrad,1);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

X = getColPoints( VHNA, 1.2, 1, 'C');
COLs=length(X);

preColMatrix=[];
elCount=1;

%create collocation matrix
for b=VHNA.el
    if messageFlag==1
       fprintf('Creating column %d of %d\n',elCount,DOFs);
    end
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
NSD=NSDlinearPhase(15);

ColMatrixNSD=zeros(size(preColMatrix));
ColRHSb=zeros(COLs,1);
tic;
for m=1:COLs
    if messageFlag==1
       fprintf('Computing row %d of %d\n',m,COLs);
    end
    for n=1:DOFs
        ColMatrixNSD(m,n)=NSD.eval(preColMatrix(m,n));
    end
    ColRHSb(m)=NSD.eval(preColRHSb(m));
end
solveTime=toc;
% fprintf('Condition number of matrix: %e',cond(ColMatrixNSD));
ColRHS=RHSa-ColRHSb;
warning('off','MATLAB:rankDeficientMatrix');
coeffs=ColMatrixNSD\ColRHS;
v_N=Projection(coeffs,VHNA);

s=linspace(0,Gamma.L,500*kwave);
semilogy(s,abs(v_N.eval(s))./kwave); xlim([Gamma.supp(1)-.1 Gamma.supp(2)+.1] ); ylim([1E-4 1E7]); 