function [ x,w ] = pathQuadV2( xi, m, S, N, freq, Rs )
%given a,b in {-\inf,0,\inf} and a weight function w(x)=exp(-x^m), for a
%function with singularities 'S', return weights and nodes

%parameter for 'being near to a singularity':
%Rs=.15;
errTol=1E-14;

%parameters for the graded quadrature:
p_max=25;   GradDelta=.15;
    
        if isempty(S)
            [x, w] = quad_gauss_exp(m, N);
        else
            if length(S)==1
                singDist = S.distFun(xi);%*freq^(1/m);
                s=S;
                %if there are two singularities... what happens here?
                if abs(singDist)<errTol && m==1 && strcmp(s.blowUpType,'log')% && 2+2==5
                    %endpoint is on singularity
                    if strcmp(s.blowUpType,'log')
                        [x, w] = quad_gengauss_loglaguerre(N);
                    else
                        error('have only coded for logarithimc singularities so far');
                    end
                elseif abs(singDist)<Rs
                    %used to be a maths-bug here for m>1
                     if m==1
                        %nearly singular, or totally singular without correct quadrature rule
                        pathSplit=sqrt(Rs^2-(singDist)^2);
                        [x1, w1]= GradedQuad( N, p_max, GradDelta );
                        %now scale the graded weights and nodes from [0 1] to [0 pathSplit]
                        x1=x1*pathSplit; w1=w1*pathSplit.*exp(-x1.^m); %HAVE JUST ADDED .^m HERE
                        %now do Gauss Laguerre on [pathSplit \infty]
                        [x2, w2] = quad_gauss_exp(m, N);
                        x2=x2+pathSplit; w2=w2*exp(-pathSplit^m);
                        %combine all weights and nodes for first SD path
                        x=[x1; x2;]; w=[w1; w2];
                     else %approximate infinite integral with finite graded integral... pretty bodgy:
                        p_max2 = 50;
                        N2 = 35;
                        [x,w] = infGradedQuad(N2, p_max2, GradDelta , m);
                     end
                else    %not even nearly singular
                    [x, w] = quad_gauss_exp(m, N);
                end
            else
                error('Have only coded for functions with single singularities');
            end
        end

end

