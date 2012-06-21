function [homolengths homodistances] = homodist(read, minlen)

% homoflags = zeros(size(read));
% % first need to identify homopolymer locations
% posL = 1;
% while posL <= length(read) - minlen +1
%         posR = posL;
%         while posR < length(read) && read(posR) == read(posR+1)
%                 posR = posR + 1;
%         end
%         homolen = posR - posL + 1;
%         if homolen >= minlen
%                 for homopos = posL:posR
%                         homoflags(homopos) = 1;
%                 end
%         end
%         posL = posR + 1;
% end

[~, homoboundaries, ~, matches] = regexp(read, '(A{3,}|T{3,}|G{3,}|C{3,})');

% % only looking for homopolymer right boundary
% homoflags(1:end-1) = (homoflags(2:end) - homoflags(1:end-1) == -1);
% homoflags(1) = 0;
% homoboundaries = find(homoflags ~= 0);
homodistances = -ones(size(read));
homolengths = homodistances;

if isempty(homoboundaries)
    return
end

for pos = 1:length(read)
        [minD, ind] = min(abs(pos - homoboundaries));
        if homoboundaries(ind) < pos  
            homodistances(pos) = minD - 1;
        else
            homodistances(pos) = minD;
        
        end
        if homodistances(pos) >= 2
            homodistances(pos) = -1;
        else
            homolengths(pos) = length(matches{ind});
        end 
end

return
