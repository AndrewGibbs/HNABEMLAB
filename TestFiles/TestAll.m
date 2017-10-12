%add all relevant paths for HNA BEM LAB:
addpath ..;
addPaths();

fileID = fopen('testData.txt','w');

%test the solutions of a simple collocation with each basis, check they
%agree:
output1=BasisTestCol( 1, 4 );
fprintf(fileID,output1);

%test the NSD quadrature
output2=QuadTestCol( 1, 6, NSDlinearPhase(15,15) );
fprintf(fileID,output2);

%display results of tests
char(output1,output2);
fclose(fileID);