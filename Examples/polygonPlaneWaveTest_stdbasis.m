clc;
clear classes;
%addPathsHNA();

%wavenumber
kwave=10;%11*2*pi;

%approximation params:
OS = 1.0; %oversampling rate
pMax=4; %polynomial degree
cL = 2;
nLayers=cL*(pMax+1)-1; %layers of mesh
hMax = 2*pi/(2*kwave);

%define the triangle
vertices =   [1    0;
              0    0;
              1/2 sqrt(3)/2;];
      
%create 'edge' object for the screen/polygon
Gamma=ConvexPolygon(vertices);

%inident plane wave
%d = [0 -1];% [-1 -1]./sqrt(2)d = [1 -1]./sqrt(2); %direction as a vector
% uinc={planeWave(kwave,[0 -1]),planeWave(kwave,[1 -1]./sqrt(2))};
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
 
%make an HNA basis on Gamma
sigmaGrad=0.15;
%V = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
V = hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);
%DOFs=length(Vstd.el);
%define the single layer 'operator' object
%A=singleLayer(kwave,Gamma);
A = combinedLayer(kwave,Gamma);

[v_h, GOA, ~, ~, T] = ColHNA(A, V, uinc, Gamma, 'oversample', OS, 'progress', 'SVDtrunc', 1E-8, 'weight','symmetry');

figure(1);
for n = 1:Gamma.numComponents
    %figure(n);
    subplot(2,2,n);
    hold on;
    s=v_h.edgeComponent(n).nodes.';
    %linspace(0, Gamma.component(n).L, 20*kwave);
    %plot(s,imag(v_h.eval(s, n)),'k', s, imag(GOA.edgeComponent(n).eval(s.')),'r');
    %plot(s,real(v_h.edgeComponent(n).eval(s)),'k', s, real(GOA.edgeComponent(n).eval(s.')),'r');
    %plot(s,imag(v_h.eval(s, n)+GOA.edgeComponent(n).eval(s.')),'r');
    %plot(s,real(GOA.edgeComponent(n).eval(s.')));
    ylim([-2*kwave 2*kwave]);
    legend('v_h');
    %hold off;
end

figure(2);
theta = linspace(0,2*pi,50*kwave);
Fv_h = FarField(Gamma, v_h, kwave, theta);
plot(theta,real(Fv_h),theta,imag(Fv_h));
%FGOA = FarField_lessSlow_stillSteady(Gamma, GOA, kwave, theta);
%plot(theta,real(Fv_N+FGOA),theta,imag(Fv_N));

%beep;