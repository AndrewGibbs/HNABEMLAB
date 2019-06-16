function [x, w] = NonOsc45(a,b,freq,N,g,singularity, width,oscs)
%alternative to NSD45 when the integrand is sufficiently non-oscillatory
    if nargin <=6
        width= b-a;
    end
    if nargin <=7
        oscs = freq*abs(b-a);
    end
    maxSingDist=.2;
    p_max=12;
    delta=.15;
    subFlag = false;
    if ~isempty(singularity)
        h = singularity.distFun;
        %else not needed
    end
    if oscs == inf
        oscs = freq*abs(b-a);
    end
    minOscsGG = 5; %Generalised Gauss cannot have arb. large number of quad points, 
                    %so cannot account for oscillations of oscillatory
                    %integral by increasing quad points
    
    Npts=max(ceil(oscs*N),N);
    
    if length(singularity)>1
        error('Can only handle one singularity per integral');
    end
    if isempty(singularity)
        [x, w0] = gradSingQuad(a, b, N, oscs, 0, delta, width);
    else
        if isempty(singularity.blowUpType)
            error('Have removed function which handles anonymous singularities');
        end
        if a < singularity.position && singularity.position < b
                [x1, w1] = NonOsc45(a,singularity.position,freq,N,g,singularity,singularity.position-a,ceil(oscs*(singularity.position-a)/width));
                [x2, w2] = NonOsc45(singularity.position,b,freq,N,g,singularity,b-singularity.position,ceil(oscs*(b-singularity.position)/width));
                x = [x1; x2];
                w0 = [w1; w2];
                subFlag = true;
        elseif strcmp(singularity.blowUpType,'log') && oscs < minOscsGG && ismember(singularity.position,[a b])
            %can use generalised Gauss quad
            if singularity.position == a
                [x, w0] = genGaussLog( Npts, a, b, width, 'L');
            elseif singularity.position == b
                [x, w0] = genGaussLog( Npts, a, b, width, 'R');
            end
        else %possible near singularity, or full singularity without GGC
            p_near = getP_near();
            if singularity.position >= b % will need to flip grading:
                [x, w0] = gradSingQuad(0, b-a, N, oscs, p_near, delta, width);
                x = b - x;
                x = flipud(x);
                w0 = flipud(w0);
            else
                [x, w0] = gradSingQuad(a, b, N, oscs, p_near, delta, width);
            end
        end
    end
    
    %only scale by oscillator on the top level of the stack
    if ~subFlag
        w = w0.*exp(1i*freq*g(x));
    else
        w = w0;
    end
    
    function p = getP_near()
        if singularity.position >= b
            p = ceil(abs(log((singularity.position-b)/(singularity.position-a))/log(delta)));
        else
            p = ceil(abs(log((a-singularity.position)/(b-singularity.position))/log(delta)));
        end
        p = min(p,p_max);
    end
    
% ----- old version: ------------------------- %
    
%     x=[];
%     if isempty(singularity)
%         [x, w0] = gradSingQuad(a, b, N, oscs, 0, delta, width);
%     elseif length(singularity)==1
%         if strcmp(singularity.blowUpType,'log') && oscs < minOscs
%             if (~isempty(singularity.position )&& oscs < minOscs) & singularity.position == a
%                 [x, w0] = genGaussLog( Npts, a, b, width, 'L');
%             elseif (~isempty(singularity.position ) && oscs < minOscs) & singularity.position == b
%                 [x, w0] = genGaussLog( Npts, a, b, width, 'R');
%             elseif ~isempty(singularity.position ) & (a < singularity.position && singularity.position < b)
%                 %split the integral and call recursively
%                 [x1, w1] = NonOsc45(a,singularity.position,freq,Npts,g,singularity);
%                 [x2, w2] = NonOsc45(singularity.position,b,freq,Npts,g,singularity);
%                 x = [x1; x2];
%                 w0 = [w1; w2];
%                 subFlag = true;
%             end
%         end
%         if isempty(x) %either not a log, or not the right type of log
%             %could be a near singularity, check for this:
%             d=width;%abs(b-a);
%             d1 = ( h(a)^2+d^2-h(b)^2 ) / (2*d);
%             d2 = ( h(b)^2+d^2-h(a)^2 ) / (2*d);
% 
%             singDist = min([h(a) h(b) sqrt(d1^2+h(a)^2)]);
%             %[singDist, nearestEndpoint] = min([abs(a-singularity.position)  abs(b-singularity.position)]);
%             if singDist<maxSingDist
%                  [x, w0] = gradSingQuad(a, b, N, oscs, p_max, delta, width);
%                 if d2<0
%                     x = b - x;
%                 end
%             else
%                 [x, w0] =  gradSingQuad(a, b, N, oscs, 0, delta, width);
%             end
%         end
%     else
%         error('Cannot do standard quadrature for multiple singularities yet');
%     end
end

