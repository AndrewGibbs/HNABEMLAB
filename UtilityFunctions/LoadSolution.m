function [vN, V] = LoadSolution( name )
%loads a previously computed solution, without the need for re-solving the
%discrete system.

    load(name);
    side=edge(vertices);
    
    if strcmp(basisType,'HNAoverlappingMesh')
        V=HNAoverlappingMesh(side, pMax, kwave, nLayers, sigmaGrad);
    elseif strcmp(basisType,'HNAsingleMesh')
        V=HNAsingleMesh(side, pMax, kwave, alphaDist, nLayers, sigmaGrad);
    elseif strcmp(basisType,'hpStandardBasis')
        V=hpStandardBasis(side, pMax, hMax, nLayers, sigmaGrad);
    else
        error('basis type not recognised');
    end
    vN=Projection(coeffs,V);

end

