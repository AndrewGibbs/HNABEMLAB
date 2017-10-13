%add all relevant paths for HNA BEM LAB:
addpath ..;
addPaths();

fileID = fopen('testData.txt','w');

%test the solutions of a simple collocation with each basis, check they
%agree:
output1=BasisTestCol( 1, 4 );
fprintf(fileID,output1);

%test the NSD quadrature for both types of HNA basis:
kwave=10; pMax=1;  nLayers=2*(pMax+1); sigmaGrad=0.15; alphaDist=2;
VHNA=HNAsingleMesh(Gamma,pMax,kwave,alphaDist, nLayers, sigmaGrad);
output2=QuadTestCol( kwave, NSDlinearPhase(15),VHNA );
VHNA=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);
output3=QuadTestCol( kwave, NSDlinearPhase(15),VHNA );
fprintf(fileID,output2);

%write a Galerkin test too, although there is no decent quadrature for this

%display results of tests
strcat(output1,output2, output3);
fclose(fileID);