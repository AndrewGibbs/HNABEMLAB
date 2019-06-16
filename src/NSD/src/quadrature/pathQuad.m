function [ x,w ] = pathQuad( a, b, m, S, N )
%given a,b in {-\inf,0,\inf} and a weight function w(x)=exp(-x^m), for a
%function with singularities 'S', return weights and nodes

%parameter for 'being near to a singularity':
Rs=.15;
errTol=1E-14;

%parameters for the graded quadrature:
p_max=10;   GradDelta=.15;

    if a==-inf && b==inf && m==2
        %away from singularities, can just do Gauss-Hermite
        [x,w] = quad_gauss_hermite(N);
    
    
    elseif a==0 && b==inf
        if isequal(S,[])
            [x, w] = quad_gauss_exp(m, N);
        else
            if length(S)==1
                s=S;
                %if there are two singularities... what happens here?
                if abs(s.position)<errTol && m==1 && strcmp(s.blowUpType,'log')
                    %endpoint is on singularity
                    if strcmp(s.blowUpType,'log')
                        [x, w] = quad_gengauss_loglaguerre(N);
                    else
                        error('have only coded for logarithimc singularities so far');
                    end
                elseif abs(s.position)<=Rs
                    %nearly singular, or totally singular without correct quadrature rule
                    pathSplit=sqrt(Rs^2-(s.position)^2);
                    [x1, w1]= GradedQuad( N, p_max, GradDelta );
                    %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
                    x1=x1*pathSplit; w1=w1*pathSplit.*exp(-x1);
                    %now do Gauss Laguerre on [pathSplit \infty]
                    [x2, w2] = quad_gauss_exp(m, N);
                    x2=x2+pathSplit; w2=w2*exp(-pathSplit^m);
                    %combine all weights and nodes for first SD path
                    x=[x1; x2;]; w=[w1; w2];
                else    %not even nearly singular
                    [x, w] = quad_gauss_exp(m, N);
                end
            else
                error('Have only coded for functions with single singularities');
            end
        end
    elseif a==0 && b<inf
        %split paths in two?
        error('Should use pathQuadFinite.m for this');
        if ~isempty(S)
            error('Havent coded Generalised Truncated Gauss-Hermite/Freud');
        end
%         [ x, w_preExp ] =  quad_gauss(N, 0, b);
%         w=w_preExp.*exp(-freq*)
        %after chatting with Daan, have decided to use standard Gauss
        %(above) rather than truncated Hermite/Frued (below)
%         if m==2
%             [ x, w ] = GaussHermiteTrunc( N, b );
%         else
%             [ x, w ] = GaussFreudTrunc( N, b, m);
%         end
        
    end

end

