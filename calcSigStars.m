function stars = calcSigStars(p)
    if p<0.001
        stars='***';
    elseif p<0.01
        stars = '**';
    elseif p<0.05
        stars = '*';
    elseif p>=0.05
        stars = 'ns';
    elseif strcmp(p,'NA')
        stars = ' ';
    end
end