function ctsandhomo2pvals(ctsandhomosfilename, cellsdir, outfilename)

% loads ctsandhomos file into memory
ctsandhomosfile = fopen(ctsandhomosfilename, 'r');
ctsandhomos = textscan(ctsandhomosfile, ...
                       '%s\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s', ...
                       'Delimiter', '\t');
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
    q_for = [ctsandhomos{8}(i) ctsandhomos{9}(i) ctsandhomos{10}(i) ctsandhomos{11}(i)];
    alphakey_for = [ctsandhomos{3}{i} '_' strrep(strrep(strrep(ctsandhomos{12}{i}, ', ', '_'), '[', ''), ']', '')];
    [p_for, alt] = polyapvalue(cts_for, alphamap(alphakey_for), ctsandhomos{3}{i});
    if isnan(p_for)
        continue
    end
    cts_rev = [ctsandhomos{16}(i) ctsandhomos{15}(i) ctsandhomos{14}(i) ctsandhomos{13}(i)];
    q_rev = [ctsandhomos{20}(i) ctsandhomos{19}(i) ctsandhomos{18}(i) ctsandhomos{17}(i)];
    alphakey_rev = [comp(ctsandhomos{3}{i}) '_' strrep(strrep(strrep(ctsandhomos{21}{i}, ', ', '_'), '[', ''), ']', '')];
    [p_rev, alt_rev] = polyapvalue(cts_rev, alphamap(alphakey_rev), comp(ctsandhomos{3}{i}));
    if isnan(p_rev)
        continue
    end
    if alt == comp(alt_rev)
        % get indices for which q score means are ref and alt
        qref_for = q_for(find(nucs == ctsandhomos{3}{i}));
        qref_rev = q_rev(find(nucs == comp(ctsandhomos{3}{i})));
        qalt_for = q_for(find(nucs == alt));
        qalt_rev = q_rev(find(nucs == comp(alt)));
        % using Fisher's method to combine p-values
        fisherXsquared = -2*(log(p_for) + log(p_rev));
        fprintf(outfile, '%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n', ...
                ctsandhomos{1}{i}, ctsandhomos{2}(i), ctsandhomos{3}{i}, alt, ...
                alphakey_for, cts_for(1), cts_for(2), cts_for(3), cts_for(4), ...
                qref_for, qalt_for, ...
                alphakey_rev, cts_rev(1), cts_rev(2), cts_rev(3), cts_rev(4), ...
                qref_rev, qalt_rev, ...
                fisherXsquared);
    end
end

fclose(outfile);
