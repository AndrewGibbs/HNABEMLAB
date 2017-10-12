%have changed single mesh and hpStandard mesh files, to include pMax on
%each mesh element. Check this doesn't mess stuff up.

clear classes;
addpath ..;
addpath ../CompGaussFiles;
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

%make an hp basis on Gamma
pMax=2; nLayers=2*pMax; sigmaGrad=0.15; hMax=2*pi/(2*kwave);
Vhp=hpStandardBasis(Gamma,pMax, hMax, nLayers, sigmaGrad);
DOFs=length(Vhp.el);

%define the single layer 'operator' object
S=SingleLayer(kwave,Gamma);

%construct Geometrical optics approximation on Gamma
GOA=GeometricalOpticsApprox(uinc,Gamma);

X = getColPoints( Vhp );
%X = ChebyshevRoots( DOFs, 'Tn', [Gamma.supp(1) Gamma.supp(2)] ).';

preColMatrix=[];
elCount=1;

%create collocation matrix
for b=Vhp.el
    Sb=S*b;
    preColMatrix=[preColMatrix Sb.col(X)];
    elCount=elCount+1;
end

%create RHSa, no integrals here, so can be computed exactly
RHS=f.eval(X);

% solve the integrals
CG=CompGaussBasic(kwave,1000,1000);

ColMatrix=zeros(size(preColMatrix));
ColRHSb=zeros(DOFs,1);
parfor n=1:DOFs
    fprintf('Computing column %d of %d\n',n,DOFs);
    for m=1:DOFs
        ColMatrix(m,n)=CG.eval(preColMatrix(m,n));
    end
end
%mix the RHSs, f(x)-A\Psi(x):
coeffs=ColMatrix\RHS;
v_N=Projection(coeffs,Vhp);
s=linspace(0,Gamma.L,100*kwave);
%plot(s,v_N.eval(s)+GOA.eval(s)); ylim([-3*kwave 3*kwave]); xlim(Gamma.supp);
plot(s,v_N.eval(s),s,GOA.eval(s)); xlim(Gamma.supp); ylim([-3*kwave 3*kwave]); 