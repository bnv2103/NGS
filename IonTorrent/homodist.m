function homodistances = homodist(read, minlen)

homoflags = zeros(size(read));
% first need to identify homopolymer locations
posL = 1;
while posL <= length(read) - minlen +1
        posR = posL;
        while posR < length(read) && read(posR) == read(posR+1)
                posR = posR + 1;
        end
        homolen = posR - posL + 1;
        if homolen >= minlen
                for homopos = posL:posR
                        homoflags(homopos) = 1;
                end
        end
        posL = posR + 1;
end

% only looking for homopolymer right boundary
homoflags(1:end-1) = (homoflags(2:end) - homoflags(1:end-1) == -1);
homoflags(1) = 0;
homoboundaries = find(homoflags ~= 0);
homodistances = length(read)*ones(size(read));
if isempty(homoboundaries)
    return
%         return homodistances
end
for pos = 1:length(read)    
        homodistances(pos) = min(abs(pos - homoboundaries));
end
            
return
