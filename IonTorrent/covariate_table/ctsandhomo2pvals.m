function ctsandhomo2pvals(ctsandhomosfilename, cellsdir, outfilename)

% loads ctsandhomos file into memory
ctsandhomosfile = fopen(ctsandhomosfilename, 'r');
ctsandhomos = textscan(ctsandhomosfile, ...
                       '%s\t%f\t%s\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f\t%f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', ...
                       'Delimiter', '\t', ...
                       'BufSize', 10000);
fclose(ctsandhomosfile);

% need to read in all alphas to one structure
alphasfiles = dir([cellsdir '/*alpha']);
alphamap = containers.Map();
for i = 1:length(alphasfiles)
    alphamap(alphasfiles(i).name(1:end-6)) = load([cellsdir '/' alphasfiles(i).name]);
end

% complement hash
comp = containers.Map({'A', 'C', 'G', 'T'}, {'T', 'G', 'C', 'A'});

% number of positions
N = length(ctsandhomos{1});

outfile = fopen(outfilename, 'w');

nucs = 'ACGT';

for i = 1:N
    cts_for = [ctsandhomos{4}(i) ctsandhomos{5}(i) ctsandhomos{6}(i) ctsandhomos{7}(i)];
    q_for = {ctsandhomos{8}{i} ctsandhomos{9}{i} ctsandhomos{10}{i} ctsandhomos{11}{i}};
    mq_for = {ctsandhomos{12}{i} ctsandhomos{13}{i} ctsandhomos{14}{i} ctsandhomos{15}{i}};
    homodat_for = eval(ctsandhomos{16}{i});
    if homodat_for(1) > 3
        homodatstrs_for = ['g3_' num2str(homodat_for(2))];
    elseif length(homodat_for) > 1
        homodatstrs_for = [num2str(homodat_for(1)) '_' num2str(homodat_for(2))];
    else
        homodatstrs_for = num2str(homodat_for);
    end
    alphakey_for = [ctsandhomos{3}{i} '_' homodatstrs_for];
    alpha_for = alphamap(alphakey_for);
    [p_for, alt] = polyapvalue(cts_for, alpha_for, ctsandhomos{3}{i});
    if isnan(p_for)
        continue
    end
    cts_rev = [ctsandhomos{20}(i) ctsandhomos{19}(i) ctsandhomos{18}(i) ctsandhomos{17}(i)];
    q_rev = {ctsandhomos{24}{i} ctsandhomos{23}{i} ctsandhomos{22}{i} ctsandhomos{21}{i}};
    mq_rev = {ctsandhomos{28}{i} ctsandhomos{27}{i} ctsandhomos{26}{i} ctsandhomos{25}{i}};
    homodat_rev = eval(ctsandhomos{29}{i});
    if homodat_rev(1) > 3
        homodatstrs_rev = ['g3_' num2str(homodat_rev(2))];
    elseif length(homodat_rev) > 1
        homodatstrs_rev = [num2str(homodat_rev(1)) '_' num2str(homodat_rev(2))];
    else
        homodatstrs_rev = num2str(homodat_rev);
    end
    alphakey_rev = [comp(ctsandhomos{3}{i}) '_' homodatstrs_rev]; 
    alpha_rev = alphamap(alphakey_rev);
    [p_rev, alt_rev] = polyapvalue(cts_rev, alpha_rev, comp(ctsandhomos{3}{i}));
    if isnan(p_rev)
        continue
    end
    if alt == comp(alt_rev)
        % get indices for which q score means are ref and alt
        refind = find(nucs == ctsandhomos{3}{i});
        qref_for = double(q_for{refind});
        mqref_for = double(mq_for{refind});
        comprefind = find(nucs == comp(ctsandhomos{3}{i}));
        qref_rev = double(q_rev{comprefind});
        mqref_rev = double(mq_rev{comprefind});
        altind = find(nucs == alt);
        qalt_for = double(q_for{altind});
        mqalt_for = double(mq_for{altind});
        compaltind = find(nucs == comp(alt));
        qalt_rev = double(q_rev{compaltind});
        mqalt_rev = double(mq_rev{compaltind});
        % skip cases with only a single nonreference, since T test gives NaN
        if cts_for(altind) < 2 || cts_rev(compaltind) < 2
            continue;
        end
        % ttest2 for q and mq on forward and reverse
        [~, q_forT] = ttest2(qref_for, qalt_for, [], 'right', 'unequal');
        [~, mq_forT] = ttest2(mqref_for, mqalt_for, [], 'right', 'unequal');
        [~, q_revT] = ttest2(qref_rev, qalt_rev, [], 'right', 'unequal');
        [~, mq_revT] = ttest2(mqref_rev, mqalt_rev, [], 'right', 'unequal');
        % fprintf(outfile, '%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        fprintf(outfile, '%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%e\t%e\t%f\t%f\t%f\t%f\t%e\t%d\t%d\t%d\t%d\t%e\t%e\t%f\t%f\t%f\t%f\t%e\n', ...
                ctsandhomos{1}{i}, ctsandhomos{2}(i), ctsandhomos{3}{i}, alt, ...
                cts_for(1), cts_for(2), cts_for(3), cts_for(4), ...
                q_forT, mq_forT, ... %log10(q_forT), log10(mq_forT), ... 
                alpha_for(1), alpha_for(2), alpha_for(3), alpha_for(4), ...
                p_for, ... %log10(p_for), ...
                cts_rev(1), cts_rev(2), cts_rev(3), cts_rev(4), ...
                q_revT, mq_revT, ... %log10(q_revT), log10(mq_revT), ...
                alpha_rev(1), alpha_rev(2), alpha_rev(3), alpha_rev(4), ...
                p_rev); %log10(p_rev));
    end
end

fclose(outfile);
