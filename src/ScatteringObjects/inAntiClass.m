function [yn] = inAntiClass(map,class,errTol)
%determines if the map between two components is inside of a given class
    if norm(map.midTranslate + class.midTranslate)<errTol && ...
            abs(map.rotate - class.rotate)<errTol && ...
            abs(map.oldLength - class.newLength)<errTol && ...
            abs(map.newLength - class.oldLength)<errTol
        yn = true;
    else
        yn = false;
    end
end