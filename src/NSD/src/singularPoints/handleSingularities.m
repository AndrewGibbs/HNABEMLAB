function [X, W] = handleSingularities(fSingularities, gSingularities, branchPoints, freq, G, N, visuals)
   if ~isempty(fSingularities)
            [fin,fon] = inpolygon(real(fSingularities),imag(fSingularities),xv,yv);
        if max(fon)
            error('f has singularity on SD path, risky business');
        end
        if max(fin)
            %got some residues to compute. Choose the radius to be small enough
            %that the integral will be non-oscillatory:
            singRad=1/freq;
            clear X_ W_;
            for s=fSingularities(fin)
                [X_, W_] = residueQuad( s,singRad, N );
                X=[X; X_]; W=[W; W_.*exp(1i*freq*G{1}(X_))];
            end
        end
   else
       X=[];    W=[];
    end
    
    badPoints=[gSingularities branchPoints];
    if ~isempty(badPoints)
        [gin,gon] = inpolygon(real(badPoints),imag(badPoints),xv,yv);
        if max(gin) || max(gon)
            error('Singularity in g(z) or branch point in region of deformation');
        end
    end
    if visuals
        plot(fSingularities,'ro','LineWidth',2);
    end
end

