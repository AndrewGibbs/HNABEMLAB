function commaCounter = nargs( f )
    %returns number of inputs to anonymous function
    fString=functions(f); 
    commaCounter=1;
    %number of input vars is number of commas plus one
   for c=1:length(fString)
      if  fString(c)==',';
          commaCounter=commaCounter+1;
      end
   end
end