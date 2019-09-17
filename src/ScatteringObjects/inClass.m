function [yn] = inClass(map,class,errTol)
%determines if the map between two components is inside of a given class
    if norm(map.midTranslate - class.midTranslate)<errTol && ...
            abs(map.rotate - class.rotate)<errTol && ...
            abs(map.oldLength - class.oldLength)<errTol && ...
            abs(map.newLength - class.newLength)<errTol
        yn = true;
    else
        yn = false;
    end
end

