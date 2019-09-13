clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
kwave=128;

%create 'screen' object ---------------------------------------------------
vertices =   [1    0;
              0    0];
Gamma=Screen(vertices);

%inident plane wave -------------------------------------------------------
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 8; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any b    asis elements
OverSample = 1.5; %choose amount to oversample by (40% here)
% construct the HNA basis (single mesh):
%VHNA = HNAsingleMesh(Gamma, pMax, kwave, throwAwayParam, nLayers, sigmaGrad, 1);
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
DOFs = length(VHNA.el); %get total #DOFs

% construct the single layer potential 'operator' ---------------------------
S=singleLayer(kwave,Gamma);

%solve (and time)
[v_N, GOA, colMatrix, colRHS, T] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress');

%plot the output
s=linspace(0,Gamma.L,min(10000,1000*kwave));
figure(1);
semilogy(s,fliplr(abs(v_N.eval(s,1))),s,fliplr(abs(v_N.evalComp(s,'+',1))),':',s,fliplr(abs(v_N.evalComp(s,'-',1))),'--','LineWidth',2);
legend('|v_N(s)|','|v_+(s)|','|v_-(s)|');
xlim([Gamma.component.supp(1)-.1 Gamma.component.supp(2)+.1] );
set(gca,'FontSize',16);
xlabel('s');
leg.FontSize = 16;
title(sprintf('k=%d',kwave),'FontSize',16);
return;%'fontsize',

%plot the far-field pattern:
figure(2);
theta = linspace(0,2*pi,500);
Fv_N = FarField(Gamma, v_N, kwave, theta);
FPsi = FarField(Gamma, GOA, kwave, theta);
plot(theta,FPsi+Fv_N);

%now plot the solution in the domain:
figure(3);
domainPlot(Gamma,uinc,GOA,v_N,kwave);
axis equal;

