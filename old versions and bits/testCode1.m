%WarningIn BEMintegral2D (line 29) : remember to add extra side for first and last side for polygons 

%basFn Warning: IF THIS EVER BECOMES A HANDLE CLASS, NEED TO CHANGE RESTRICTION BELOW

%for polygons, will need to change:
    %edge.m trace/normal functions
    %make sure normals point outwards

clear classes;
%wavenumber
kwave=5;

%define the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
      
%create 'edge' object for the screen
Gamma=edge(vertices);

%make an HNA basis on Gamma
pMax=2; nLayers=1*pMax; sigmaGrad=0.15;
V=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);

%define the single layer 'operator' object
S_k=2*SingleLayer(kwave);

% %initialise Galerkin matrix
% GalerkinMatrix=zeros(V.numEls);

for m=1:V.numEls
    for n=1:V.numEls
        GalerkinMatrix{n,m}=BEMintegral2D(S_k, V.el(m), V.el(n));
        %(note there was a bizzare error I couldn't understand when I tried
        %to write this as a matrix, not a structure)
    end
end