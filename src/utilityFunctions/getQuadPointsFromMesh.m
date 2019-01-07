function [s,w] = getQuadPointsFromMesh(meshIn,pointsPerWavelength,k,distType)
 % distributed by distType = 'U' for uniform, 'G' for Gauss
        s = []; w=[];
   for n = 1:length(meshIn.el)
                
                if strcmp(distType,'U')
                    s_ = linspace(meshIn.points(n),meshIn.points(n+1),N).';
                    if n<length(meshIn.el)
                        s_ = s_(1:(end-1));
                    end
                elseif strcmp(distType,'G')
                    %[s_,w_] =gauss_quad(self.points(n),self.points(n+1),N);
                     [s_,w_] = gauss_quad_wave_split2(meshIn.points(n), meshIn.points(n+1),...
                                pointsPerWavelength, k,  meshIn.el(n).width );
                end
                s = [s; s_];
                w = [w; w_];
                clear s_ w_;
   end
end

