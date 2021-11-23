clear all;
%wavenumber
kwave=4*pi;%11*2*pi;

%approximation params:
OS = 1.5; %oversampling rate
pMax=4; %polynomial degree
cL = 2;
nLayers=cL*(pMax+1)-1; %layers of mesh
% hMax = 2*pi/(2*kwave);

%define the square
%vertices = [1 0; 0 -1; -1 0; 0 1];
vertices = sqrt(2)*[-.5 -.5; -.5 .5; .5 .5; .5 -.5];
      
%create 'edge' object for the screen/polygon
Gamma=ConvexPolygon(vertices);

%inident plane wave
d = [1 0];
uinc=planeWave(kwave,d);
 
%make an HNA basis on Gamma
sigmaGrad=0.15;
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
hmax=1/kwave;
% VHNA = hpStandardBasis(Gamma, pMax, hmax, nLayers, sigmaGrad);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
%A=singleLayer(kwave,Gamma);
A = combinedLayer(kwave,Gamma);

[v_HNA, GOA, colMatrix, colRHS] = ColHNA(A, VHNA, uinc, Gamma, 'oversample', OS, 'progress', 'SVDtrunc', 1E-8, 'weight');

figure(2);
domainPlot(Gamma,uinc,GOA,v_HNA,kwave,10000);


figure(3);

theta = linspace(0,2*pi,50*kwave);
Fv_N = FarField(Gamma, v_HNA, kwave, theta);
FGOA = FarField(Gamma, GOA, kwave, theta);
plot(theta,real(Fv_N+FGOA),theta,imag(Fv_N+FGOA),theta,abs(Fv_N+FGOA));
% plot(theta,real(Fv_N),theta,imag(Fv_N),theta,abs(Fv_N));
