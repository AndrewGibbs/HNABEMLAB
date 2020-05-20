function I = RHSquad(Operator, GOA, edge_index, Xcol, Nquad)
%computes oversampled collocation projection using HNA basis/frame
    if isa(GOA.uinc,'planeWave')
        if edge_index == Xcol.side
            I = colEvalLinear(Operator, GOA, edge_index, Xcol, Nquad, []);
        else
            I = bruteForce(Operator, GOA, edge_index, Xcol);
        end
    else
        error('Havent coded for other wave types yet');
    end 
end