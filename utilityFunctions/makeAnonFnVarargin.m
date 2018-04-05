function phaseText= makeAnonFnVarargin( LHS, FnName, dim )

       phaseText=sprintf(' @(varargin) %s(varargin{1}',FnName);
       for n=2:dim
           phaseText=strcat(phaseText, sprintf(',varargin{%d}',n));
       end
       phaseText=strcat(phaseText, ');');
       %eval(phaseText);
end

