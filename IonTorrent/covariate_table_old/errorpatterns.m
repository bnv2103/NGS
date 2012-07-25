function errorpatterns(fasta, sam)

%[reffile,refpath] = uigetfile({'*.fasta;*.fa'}, ...
%                                          'select reference fasta file');
%[fahdr, refseq] = fastaread([refpath, reffile]);
disp('loading fasta')
[fahdr, refseq] = fastaread(fasta, 'TrimHeaders', true, 'Blockread', [1 24]);
%refseq = FASTAstruct.Sequence;
%[alnfile,alnpath] = uigetfile({'*.sam'}, ...
%                                          'select alignment sam file');
%[SAMstruct, SAMhdr] = samread([alnpath, alnfile]);
disp('loading sam')
[SAMstruct, ~] = samread(sam);
unmapped_flags = bitand(bitshift(cell2mat({SAMstruct.Flag}), -2), 1);
SAMstruct = SAMstruct(unmapped_flags == 0); % exclude unmapped reads

% % debug code
% whos
% size(refseq)
% %

minhomopolymerlength = 3;

reads = {SAMstruct.Sequence};
positions = cell2mat({SAMstruct.Position});
qualities = {SAMstruct.Quality};
cigars = {SAMstruct.CigarString};
reverse_flags = double(bitand(bitshift(cell2mat({SAMstruct.Flag}), ...
                                       -4), ...
                              1)...
                      );
chromosomes = {SAMstruct.ReferenceName};
                  
clear FilterIndex fahdr unmapped_flags SAMstruct SAMhdr
% this shows deletions as hyphens and omits insertsions. Later if we want
% to include indels we can set GapsInRef to true to show insertions
% assuming hyphens in refernce at insert position (different for each read)
% note that we also are losing pad/clip sequence
nreads = numel(reads);

%f = fopen([alnpath, alnfile, '.errorpatterns'], 'w');
f = fopen([sam, '.errorpatterns'], 'w');

chrs = {'chr1'
        'chr2'
        'chr3'
        'chr4'
        'chr5'
        'chr6'
        'chr7'
        'chr8'
        'chr9'
        'chr10'
        'chr11'
        'chr12'
        'chr13'
        'chr14'
        'chr15'
        'chr16'
        'chr17'
        'chr18'
        'chr19'
        'chr20'
        'chr21'
        'chr22'
        'chrX'
        'chrY'};

for i = 1:nreads
    chr = chromosomes{i};
    if ~ismember(chr, chrs)
        continue
    end
    
    
    
    aligned_read = cigar2align(reads(i), cigars(i));
    pos = positions(i) + 1:length(aligned_read) - 1;
    aligned_quality = cigar2align(qualities(i), cigars(i));
    aligned_ref  = refseq{strcmp(fahdr, chr)}(positions(i):positions(i) ...
                          + length(aligned_read) - 1);
    % reverse complement
    if reverse_flags(i)
        aligned_read = seqrcomplement(aligned_read);
        aligned_ref = seqrcomplement(aligned_ref);
        aligned_quality = fliplr(aligned_quality);
    end
    % ignore deletions from ref (insertions were removed by cigar2align)
    dels = aligned_read ~= '-';
    aligned_read = aligned_read(dels);
    aligned_ref = aligned_ref(dels);
    aligned_quality = aligned_quality(dels);
    pos = pos(dels);
    
    % compute the error pattern variables
    q = double(aligned_quality);
    x = 1:length(aligned_read);
    [hs, hd] = homodist(aligned_read, minhomopolymerlength);
    GC = (length(strfind(aligned_read, 'G', 1)) + length(strfind(aligned_read, 'G', 2))) / length(aligned_read); 

    for n = 1:x(end)
        fprintf(f, '%s\t%u\t%c\t%c\t%u\t%u\t%u\t%u\t%u\t%u\n', ...
                                               chr, ...
                                               pos(n), ...
                                               aligned_ref(n), ...
                                               aligned_read(n), ...
                                               reverse_flags(i), ...
                                               q(n), ...
                                               x(n), ...
                                               hs(n), ...
                                               hd(n), ...
                                               GC);
    end
    
    clc
    fprintf('%u of %u reads processed', i, nreads);
end

fclose(f);
