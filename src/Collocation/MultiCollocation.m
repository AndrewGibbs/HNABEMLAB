function self=MultiCollocation(copyCol, Vbasis, self)
    %copyBasis must be a function handle of the form:
    % @(x) aBasis(x, ...) where x is the obstacle, and aBasis is
    % some subclass of Basis
    %self.obstacle=Vbasis.obstacle;

    onSide = [];
    cumNumPtsPrev = 1;
    for n = 1:length(Vbasis.edgeBasis)
         self.edgeCol{n} = copyCol(Vbasis.edgeBasis{n});
         cumNumPtsNext = length(self.edgeCol{n}.pt)+cumNumPtsPrev-1;
         self.pt(cumNumPtsPrev:cumNumPtsNext) = self.edgeCol{n}.pt;
         self.onEdge(cumNumPtsPrev:cumNumPtsNext) = n;
         cumNumPtsPrev = cumNumPtsNext + 1;
    end
     
    self.numPoints = cumNumPtsNext;
    self.ruleType = self.edgeCol{1}.ruleType;
end