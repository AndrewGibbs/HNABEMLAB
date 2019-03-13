
function self=MultiBasis(copyBasis, obstacle, self)
    %copyBasis must be a function handle of the form:
    % @(x) aBasis(x, ...) where x is the obstacle, and aBasis is
    % some subclass of Basis
    self.obstacle=obstacle;

    indexEnd = 0;
    self.meshDOFs = [];
    for n=1:obstacle.numComponents
        indexStart = indexEnd + 1;
        self.edgeBasis{n} = copyBasis(obstacle.component{n});
        indexEnd = indexStart -1 + length(self.edgeBasis{n}.el);
        self.mesh{n} = self.edgeBasis{n}.mesh;

        self.plusCoefs    = [self.plusCoefs    indexStart+self.edgeBasis{n}.plusCoefs];
        self.minusCoefs   = [self.minusCoefs   indexStart+self.edgeBasis{n}.minusCoefs];
        self.nonOscCoeffs = [self.nonOscCoeffs indexStart+self.edgeBasis{n}.nonOscCoeffs];
        self.el           = [self.el           self.edgeBasis{n}.el];
        self.meshDOFs     = [self.meshDOFs     self.edgeBasis{n}.meshDOFs];
        self.elSide(indexStart:indexEnd) = n;
    end

end

