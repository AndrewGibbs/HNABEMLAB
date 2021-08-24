% clc;
% clear classes;
%addPathsHNA();

%wavenumber
kwave=10;%11*2*pi;

%approximation params:
OS = 1.5; %oversampling rate
pMax=4; %polynomiail degree
cL = 2;
nLayers=cL*(pMax+1)-1; %layers of mesh
% hMax = 2*pi/(2*kwave);

%define the triangle
vertices =   [1    0;
              0    0;
              1/2 sqrt(3)/2;];
      
%create 'edge' object for the screen/polygon
Gamma=ConvexPolygon(vertices);

%inident plane wave
d = [1 -1]./sqrt(2);% [-1 -1]./sqrt(2)
uinc=planeWave(kwave,d);
 
%make an HNA basis on Gamma
sigmaGrad=0.15;
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
%VHNA = hpStandardBasis(Gamma, pMax, hMax, nLayers, sigmaGrad);
DOFs=length(VHNA.el);
%define the single layer 'operator' object
%A=singleLayer(kwave,Gamma);
A = combinedLayer(kwave,Gamma);

[v_HNA, GOA, colMatrix, colRHS] = ColHNA(A, VHNA, uinc, Gamma, 'oversample', OS, 'progress', 'SVDtrunc', 1E-8, 'weight');

% for n = 1:Gamma.numComponents
%     s=linspace(0, Gamma.component(n).L, 20*kwave);
%     figure(n);
%     semilogy(s,abs(v_HNA.eval(s,n))); %ylim([1E-4 1E7]); 
% end
%figure(n+1);
for n = 1:Gamma.numComponents
    %figure(n);
    subplot(2,2,n);
    hold on;
    s=linspace(0, Gamma.component(n).L, 20*kwave);
        plot(s,real(v_HNA.evalComp(s, '+', n)));
        plot(s,real(v_HNA.evalComp(s, '-', n)));
        plot(s,real(v_HNA.eval(s, n)),'k');
        plot(s,real(GOA.edgeComponent(n).eval(s.')));
    ylim([-2*kwave 2*kwave]);
    legend('v_+','v_-','v_+e_++v_-e_-','GOA');
    hold off;
end

figure(3);
theta = linspace(0,2*pi,50*kwave);
Fv_N = FarField(Gamma, v_HNA, kwave, theta);
FGOA = FarField(Gamma, GOA, kwave, theta);
plot(theta,real(Fv_N+FGOA),theta,imag(Fv_N+FGOA));
%plot(theta,real(Fv_N),theta,imag(Fv_N));

beep;
% for n = 1:Gamma.numComponents
%     figure(1);
%     subplot(2,2,n);
%     hold on;
%     s=linspace(0, Gamma.component(n).L, 20*kwave);
%         plot(s,real(v_HNA.eval(s, n)),'k');
%         plot(s,real(GOA.edgeComponent(n).eval(s.')));
%     ylim([-20*kwave 20*kwave]);
%     legend('v_N','GOA');
%     hold off;
% end