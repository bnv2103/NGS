[reffile,refpath,FilterIndex] = uigetfile({'*.fasta;*.fa'}, ...
                                          'select reference fasta file');
[fahdr, refseq] = fastaread([refpath, reffile]);
[alnfile,alnpath,FilterIndex] = uigetfile({'*.sam'}, ...
                                          'select alignment sam file');
[SAMstruct, SAMhdr] = samread([alnpath, alnfile]);
unmapped_flags = bitand(bitshift(cell2mat({SAMstruct.Flag}), -2), 1);
SAMstruct = SAMstruct(unmapped_flags == 0); % exclude unmapped reads

minhomopolymerlength = 4;

reads = {SAMstruct.Sequence};
positions = cell2mat({SAMstruct.Position});
qualities = {SAMstruct.Quality};
cigars = {SAMstruct.CigarString};
reverse_flags = double(bitand(bitshift(cell2mat({SAMstruct.Flag}), ...
                                       -4), ...
                              1)...
                      );
clear FilterIndex fahdr unmapped_flags SAMstruct SAMhdr
% this shows deletions as hyphens and omits insertsions. Later if we want
% to include indels we can set GapsInRef to true to show insertions
% assuming hyphens in refernce at insert position (different for each read)
% note that we also are losing pad/clip sequence
nreads = numel(reads);

f = fopen([alnpath, alnfile, '.errorpatterns'], 'w');

for i = 1:nreads
    aligned_read = cigar2align(reads(i), cigars(i));
    aligned_quality = cigar2align(qualities(i), cigars(i));
    aligned_ref  = refseq(positions(i):positions(i) ...
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
    
    % compute the error pattern variables
    q = double(aligned_quality);
    x = 1:length(aligned_read);
    h = homodist(aligned_read, minhomopolymerlength);
    
    for n = 1:x(end)
        fprintf(f, '%c\t%c\t%u\t%u\t%u\t%u\n', aligned_ref(n), ...
                                               aligned_read(n), ...
                                               reverse_flags(i), ...
                                               q(n), ...
                                               x(n), ...
                                               h(n));
    end
    
    clc
    fprintf('%u of %u reads processed', i, nreads);
end

fclose(f);