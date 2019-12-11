clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
%wavelengthsPerComponent = 100;
%scatter by Cantor set of order
CantOrder = 2;
L = 1;
% now choose the wavenumber based on above info --------------------------

%kwave = ceil(2*pi*wavelengthsPerComponent*3^CantOrder);

kwave = 128;

fprintf('wavenumber chosen to be k=%d',kwave);

%create 'screen' object ---------------------------------------------------
vertices =   [L    0;
              0    0];
          
segs = Cantor(CantOrder,L);

Gamma=MultiScreen(vertices,segs);

%inident plane wave -------------------------------------------------------
%d = [-1 0]; %direction as a vector
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 4 ; %polynomial degree
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
[v_N, GOA, colMatrix, colRHS] = ColHNA(S, VHNA, uinc, Gamma,'oversample', OverSample, 'progress','quadpoints',50);
T = toc;

disp('Plotting output');

%now plot the solution in the domain:
figure(1);
domainPlot(Gamma,uinc,GOA,v_N,kwave);
axis equal;
ylim([-.2 .2]);
title(sprintf('k = %d',kwave));