clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
kwave=250;

%create 'screen' object ---------------------------------------------------
vertices =   [1    0;
              0    0];
          
segs = [0 .45 .5 1];

Gamma=MultiScreen(vertices,segs);

%inident plane wave -------------------------------------------------------
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 16 ; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.5; %choose amount to oversample by (50% here)
% construct the HNA basis (single mesh):
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
%VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
DOFs = length(VHNA.el); %get total #DOFs

% construct the single layer potential 'operator' ---------------------------
S=singleLayer(kwave,Gamma);

%solve (and time)
tic;
[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress');
T = toc;


disp('Plotting output');

%plot the far-field pattern:
figure(1);
theta = linspace(0,2*pi,5000);
Fv_N = FarField(Gamma, v_N, kwave, theta);
FPsi = FarField(Gamma, GOA, kwave, theta);
plot(theta,(FPsi(:,1)+Fv_N(:,1)));
xlim([0 2*pi]);
xlabel('\theta');

%now plot the solution in the domain:
figure(2);
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

%now compute the complimentary sound-hard aperature problem, by Babinet's
%principle:
figure(3);
yU = linspace(0,ymax,pixelRate*(ymax));
yL = linspace(ymin,0,pixelRate*(-ymin));
yL = yL(1:(end-1));
SvU = singleLayerDomain(Gamma, v_N, kwave, x, yU);
SPsiU = singleLayerDomain(Gamma, GOA, kwave, x, yU);
SvL = singleLayerDomain(Gamma, v_N, kwave, x, yL);
SPsiL = singleLayerDomain(Gamma, GOA, kwave, x, yL);

uincR = uinc.getReflection(Gamma);
[X1U, X2U] = meshgrid(x,yU);
u_NU = uinc.eval(X1U,X2U) - SvU.' - SPsiU.';

DirZone = .01/kwave;
B_U = uincR.eval(X1U,X2U) + forceDirichlet(x,yU,u_NU,DirZone,Gamma);

[X1L, X2L] = meshgrid(x,yL);
u_NL = uinc.eval(X1L,X2L) - SvL.' - SPsiL.';
B_L = uinc.eval(X1L,X2L) - forceDirichlet(x,yL,u_NL,DirZone,Gamma);

babz = [B_L; B_U;];
%babz = [forceDirichlet(x,yL,B_L,DirZone,Gamma); forceDirichlet(x,yU,B_U,DirZone,Gamma);];
imagesc(x,fliplr([yL yU]),flipud(real(babz)));

shading interp;
set(gca,'YDir','normal');
hold on;
Gamma.drawCompliment(xmin,xmax);