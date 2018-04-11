function p_i = pMaxChoose( index, pMax, nStar )
    %based on HNA 2013 paper remark 5.3
    if index<nStar
        p_i=ceil(pMax*(index-1)/nStar);
    else
        p_i=pMax;
    end
end

