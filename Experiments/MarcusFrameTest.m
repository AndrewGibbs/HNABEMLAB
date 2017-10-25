%-------------------------------------------------------------------------%
kwave=1000;                                         %wavenumber of problem
pMax=10;                                %polynomial degree of HnA space
nLayers=2*(pMax+1); sigmaGrad=0.15;         %grading params for corners
extraPoints=1;                       %extraPonts per DOF on each element
%-------------------------------------------------------------------------%

%create 'edge' object for the screen
vertices=[0 0       %first vertex
          1 0];     %second vertex
Gamma=edge(vertices);

VHNA=HNAoverlappingMesh(Gamma,pMax,kwave, nLayers, sigmaGrad);
basEls=length(VHNA.el);

X = getColPoints( VHNA, extraPoints, 1);
   
for m=1:basEls
    for n=1:length(X)
        MarcusMatrix(m,n)=VHNA.el(m).eval(X(n));
    end
end

MarcusSingularValues=svd(MarcusMatrix);

semilogy(MarcusSingularValues,'x')