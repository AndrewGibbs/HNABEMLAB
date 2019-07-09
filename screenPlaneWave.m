clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
kwave=50;

%create 'screen' object ---------------------------------------------------
vertices =   [1    0;
              0    0];
Gamma=Screen(vertices);

%inident plane wave -------------------------------------------------------
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 14; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.5; %choose amount to oversample by (40% here)
% construct the HNA basis (single mesh):
%VHNA = HNAsingleMesh(Gamma, pMax, kwave, throwAwayParam, nLayers, sigmaGrad, 1);
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
DOFs = length(VHNA.el); %get total #DOFs

% construct the single layer potential 'operator' ---------------------------
S=singleLayer(kwave,Gamma);

%solve (and time)
tic;
[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress');
T = toc;

%plot the output
s=linspace(0,Gamma.L,min(1000,10000*kwave));
figure(1);
semilogy(s,abs(v_N.eval(s,1))); ylim([1E-2 1E3]);
xlim([Gamma.component.supp(1)-.1 Gamma.component.supp(2)+.1] );

%plot the far-field pattern:
figure(2);
theta = linspace(0,2*pi,500);
Fv_N = FarField(Gamma, v_N, kwave, theta);
FPsi = FarField(Gamma, GOA, kwave, theta);
plot(theta,FPsi+Fv_N);

%now plot the solution in the domain:
figure(3);
totalPixels = 1000000;
xmin = -.3; xmax = 1.3;
ymin = -.2; ymax = .2;
imageArea = (xmax-xmin)*(ymax-ymin);
pixelRate = ceil(sqrt(totalPixels/imageArea));
figure(2);

y = linspace(ymin,ymax,pixelRate*(ymax-ymin));
x = linspace(xmin,xmax,pixelRate*(xmax-xmin));
Sv = singleLayerDomain(Gamma, v_N, kwave, x, y);
SPsi = singleLayerDomain(Gamma, GOA, kwave, x, y);
[X1, X2] = meshgrid(x,y);
u_N = uinc.eval(X1,X2) - Sv.' - SPsi.';
imagesc(x,fliplr(y),flipud(real(u_N)));
shading interp;
set(gca,'YDir','normal');
hold on;
Gamma.draw;

