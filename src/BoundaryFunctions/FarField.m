function F = FarField( self  )
%create the far-field integral, but do not evaluate.
    boundary{1}=[0 2*pi];
    boundary{2}=self.domain;
    
    F=BoundaryIntegral(kernel, self, boundary);

end

