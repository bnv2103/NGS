

samfile = uigetfile('*.sam');
[SAMstruct, ~] = samread(samfile);
unmapped_flags = bitand(bitshift(cell2mat({SAMstruct.Flag}), -2), 1);
SAMstruct = SAMstruct(unmapped_flags == 0); % exclude unmapped reads


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

nreads = numel(reads);

f = fopen([samfile, '.alignedMATLAB'], 'w');

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
    if ~ismember(chr, chrs) || positions(i) < 1 || strcmp(reads(i), '*') || strcmp(qualities(i), '*')
        continue
    end
    
    aligned_read = cigar2align(reads(i), cigars(i));
    aligned_quality = cigar2align(qualities(i), cigars(i));  
    fprintf(f, '%s\t%s\n', aligned_read, aligned_quality);
    
end
    
fclose(f);
    
