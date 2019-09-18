clc;
clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
wavelengthsPerComponent = 100;
%scatter by Cantor set of order
CantOrder = 5;

% now choose the wavenumber based on above info --------------------------

kwave = ceil(2*pi*wavelengthsPerComponent*3^CantOrder);
fprintf('wavenumber chosen to be k=%d',kwave);

%create 'screen' object ---------------------------------------------------
vertices =   [1    0;
              0    0];
          
segs = Cantor(CantOrder,1);

Gamma=MultiScreen(vertices,segs);

%inident plane wave -------------------------------------------------------
d = [1 -1]./sqrt(2); %direction as a vector
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 8 ; %polynomial degree
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

theta = linspace(0,2*pi,20*kwave);
if kwave<4000
    Fv_N = FarField(Gamma, v_N, kwave, theta);
    FPsi = FarField(Gamma, GOA, kwave, theta);
else
    Fv_N = FarField_slowNsteady(Gamma, v_N, kwave, theta);
    FPsi = FarField_slowNsteady(Gamma, GOA, kwave, theta);
end
FF = Fv_N + FPsi;
polarplot(theta,log10(abs(flipud(FF))));
radMin = -1.5;
hold on;
rlim([radMin,max(log10(abs(FF)))]);
segsShitft = abs(radMin)*(segs-.5)+radMin;
for n = 1:length(segs)/2
    polarplot([0; 0],[segsShitft(2*n-1); segsShitft(2*n)],'k','LineWidth',2);
end
hold on;

return;
%now plot the solution in the domain:
figure(2);
domainPlot(Gamma,uinc,GOA,v_N,kwave);
axis equal;

%now compute the complimentary sound-hard aperature problem, by Babinet's
%principle:
figure(3);
BabinetComplementPlot(Gamma,uinc,GOA,v_N,kwave);
axis equal;