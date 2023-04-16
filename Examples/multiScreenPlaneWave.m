clear classes;
% run addPathsHNA() to add necessary search paths
%wavenumber
%wavelengthsPerComponent = 100;
%scatter by Cantor set of order
%CantOrder = 1;

% now choose the wavenumber based on above info --------------------------

%kwave = ceil(2*pi*wavelengthsPerComponent*3^CantOrder);

kwave = 30;

fprintf('wavenumber chosen to be k=%d',kwave);

%create 'screen' object ---------------------------------------------------
vertices =   [1    0;
              0    0];
          
CantOrder = 5;
segs = Cantor(CantOrder,1/3);
%segs = 10*pi-[0 2*pi 21*pi/10 5*pi/2 14*pi/5 7*pi/2 4*pi 6*pi 61*pi/10 10*pi];
%segs = [0 0.1 0.3 0.5 0.7 1];

Gamma=MultiScreen(vertices,segs);

%inident plane wave -------------------------------------------------------
d = [0.5000 -0.8660]; %direction as a vector
d = d/norm(d);
uinc=planeWave(kwave,d);
    
%make an HNA basis on Gamma -----------------------------------------------
pMax = 4 ; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.5; %choose amount to oversample by (50% here)
% construct the HNA basis (overlapping mesh):
VHNA = HNAoverlappingMesh(Gamma, pMax, kwave, nLayers, sigmaGrad);
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

theta = linspace(0,2*pi,50*kwave);%5000 seems apt
theta = theta(1:(end-1)); %delete final point here
Fv_N = FarField(Gamma, v_N, kwave, theta);
FPsi = FarField(Gamma, GOA, kwave, theta);
FF = Fv_N + FPsi;
radMin = -2;
FFlogMin = -1;
FFlog = log10(abs((FF)));
logAbsFFhandle = @(theta_) log10(abs(FarField(Gamma, v_N, kwave, theta_) + FarField(Gamma, GOA, kwave, theta_)));
%now adjust (theta,FF) to include trophs in log plot
[theta, FFlog] = getFFtroughs(theta, FFlog, logAbsFFhandle);
FFlog = flipud(FFlog);
FFlog(FFlog<FFlogMin) = FFlogMin;
polarplot(theta,FFlog);
% hold on;
% rlim([radMin,3.5]);
% segsShitft = abs(radMin)*(segs-.5)+radMin;
% for n = 1:length(segs)/2
%     polarplot([0; 0],[segsShitft(2*n-1); segsShitft(2*n)],'k','LineWidth',2);
% end
% 
return;
%now plot the solution in the domain:
figure(2);
domainPlot(Gamma,uinc,GOA,v_N,kwave,10000);
% axis equal;

%now compute the complimentary sound-hard aperature problem, by Babinet's
%principle:
% figure(3);
% BabinetComplementPlot(Gamma,uinc,GOA,v_N,kwave);
% axis equal;