function  SaveSolution( V,coeffs,name,vertices,uinc )
    if nargin>=5
       uincType=class(uinc);
       if isequal(uincType,'planeWave')
           d=uinc.d;
       else
           warning('need to code this still');
       end
    end
    basisType=class(V);
    pMax=V.pMax;
    if strcmp(class(V),'hpStandardBasis')
        hMax=V.hMax;
    elseif strcmp(class(V),'HNAsingleMesh')
        alphaDist=V.alphaDist;
    end
    nLayers=V.nLayers;
    sigmaGrad=V.sigmaGrad;
    kwave=uinc.kwave;
    %clear the objects, these don't save well
    clear V uinc;
    %save what's left, which is enough to reconstruct the basis
    save(name);
end

