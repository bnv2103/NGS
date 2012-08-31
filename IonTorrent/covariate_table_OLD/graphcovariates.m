f = fopen('sample.sam.errorpatternsNOGERM', 'r');
C = textscan(f, '%s\t%s\t%c\t%c\t%u\t%u\t%u\t%u\t%d\t%d\t%f\n');
fclose(f);

% ref = C{3};
% read = C{4};
% len = C{5};
% reverse = logical(C{6});
% quality = C{7};
% position = C{8};
% hs = C{9};
% hd = C{10};
% GC = C{11};

N = length(C{5});

match = (C{3} == C{4});

ref = C{3};
As = (C{4} == 'A');
Cs = (C{4} == 'C');
Gs = (C{4} == 'G');
Ts = (C{4} == 'T');

spc = 30;
% l = linspace(double(min(C{5})), double(max(C{5})), spc);
% p = linspace(1, double(max(C{5})), spc);
q = double(min(C{7})):double(max(C{7}));
% gc = linspace(double(min(C{11})), double(max(C{11})), spc);
hs = 1:double(max(C{9}));
hd = [-1 0 1];
% 
% % disp('reverse position (from 3'' end)')
% % l_cts = [hist(C{5}, l); hist(C{5}(~logical(match)), l)];
% % p_cts = [hist(C{5}-C{8}, p); hist(C{5}(~logical(match))-C{8}(~logical(match)), p)];
% q_cts = [histc(C{7}, q); histc(C{7}(~logical(match)), q)];
% % gc_cts = [hist(C{11}, gc); hist(C{11}(~logical(match)), gc)];
% hs_cts_0 = [histc(C{9}(logical(C{10} == 0)), hs); histc(C{9}(~logical(match) & logical(C{10} == 0)), hs)];
% hs_cts_1 = [histc(C{9}(logical(C{10} == 1)), hs); histc(C{9}(~logical(match) & logical(C{10} == 1)), hs)];
% 
% thresh = 0;%100000;
% 
% figure(1)
% % subplot(2,2,1)
% % plot(l(l_cts(1,:)>thresh), l_cts(2,l_cts(1,:)>thresh)./l_cts(1,l_cts(1,:)>thresh))
% % xlabel('read length')
% % 
% % subplot(2,2,2)
% % plot(p(p_cts(1,:)>thresh), p_cts(2,p_cts(1,:)>thresh)./p_cts(1,p_cts(1,:)>thresh))
% % xlabel('base position in reaad')
% 
% % subplot(2,2,3)
% % plot(10.^(-(q(q_cts(1,:)>thresh) - 33)/10), q_cts(2,q_cts(1,:)>thresh)./q_cts(1,q_cts(1,:)>thresh), '--ko')
% scatter(q(q_cts(1,:)>thresh) - 33, ...
%         -10*log10(q_cts(2,q_cts(1,:)>thresh)./q_cts(1,q_cts(1,:)>thresh)), ...
%         .000035*(q_cts(1,q_cts(1,:)>thresh)), 'filled')
% xlabel('base quality')
% ylabel('phred scaled mismatch ratio')
% axis([0 40 0 40])
% axis square
% hold on
% plot([0 40], [0 40], 'k')
% mean_mismatch = -10*log10(sum(~match)/length(match));
% scatter([mean_mismatch], [mean_mismatch], 100, 'rs', 'filled')
% hold off
% % subplot(2,2,4)
% % plot(gc(gc_cts(1,:)>thresh), gc_cts(2,gc_cts(1,:)>thresh)./gc_cts(1,gc_cts(1,:)>thresh))
% % xlabel('read GC content')
% 
% figure(2)
% scatter(hs(hs_cts_0(1,:)>thresh), -10*log10(hs_cts_0(2,hs_cts_0(1,:)>thresh)./hs_cts_0(1,hs_cts_0(1,:)>thresh)), ...
%     .001*(hs_cts_0(1,hs_cts_0(1,:)>thresh)), 'filled', 'r')
% hold on
% scatter(hs(hs_cts_1(1,:)>thresh), -10*log10(hs_cts_1(2,hs_cts_1(1,:)>thresh)./hs_cts_1(1,hs_cts_1(1,:)>thresh)), ...
%          .001*(hs_cts_1(1,hs_cts_1(1,:)>thresh)), 'filled', 'b')
% hold off
% xlabel('homopolymer length')
% legend('at homopolymer 3'' boundary', '1 basepair from boundary')
% ylim([0 40])

matchpoints = sortrows([int32(C{7}(logical(match))) C{9}(logical(match),:) C{10}(logical(match),:)]);
nmatchpoints = size(matchpoints, 1);
mismatchpoints = sortrows([int32(C{7}(~logical(match))) C{9}(~logical(match),:) C{10}(~logical(match),:)]);
nmismatchpoints = size(mismatchpoints, 1);

% cells for all points, match and mismatch
cells = sortrows(unique([matchpoints; mismatchpoints], 'rows'));
ncells = size(cells, 1);

match_cts = zeros(ncells, 1);
mismatch_cts = zeros(ncells, 1);

m = 1;
mm = 1;
for c = 1:ncells
%     for m = 1:nmatchpoints
    while m <= nmatchpoints && all(cells(c, :) == matchpoints(m, :))
        match_cts(c) = match_cts(c) + 1;
        m = m + 1;
    end
    while mm <= nmismatchpoints && all(cells(c, :) == mismatchpoints(mm, :))
        mismatch_cts(c) = mismatch_cts(c) + 1;
        mm = mm + 1;
    end
    clc; disp(100*c/ncells);
end

% figure(3)
% edges = linspace(0, 7, 35);
% bar(edges, [histc(log10(match_cts), edges) histc(log10(mismatch_cts), edges)])

dlmwrite('sample.sam.covariatetable', [cells match_cts mismatch_cts], '\t')
