[FileName,PathName,FilterIndex] = uigetfile('*.sam','select reads SAM file');
disp('loading reads from SAM file...')
reads = samread([PathName,FileName]);
ref_string = unique({reads.ReferenceName});
ref_string = ref_string(find(~strcmp(ref_string,'*')));
assert(length(ref_string) == 1);
reads = reads(strcmp(ref_string, {reads.ReferenceName}));
reads = BioMap(reads);
starts = getStart(reads);
reads = cigar2align(reads.Sequence, reads.Signature);

nreads = size(reads,1);

[FileName,PathName,FilterIndex] = uigetfile('*.fasta','select reference FASTA file');
disp('loading reference from fasta file...')
reference = fastaread([PathName,FileName]);
reference = reference.Sequence;

clear ref_string FileName PathName FilterIndex

disp('computing mismatches...')
max_read_len = 0;
for i = 1:nreads
    read = deblank(reads(i,:));
    len = length(read);s
    if(len > max_read_len) max_read_len = len; end
end

mismatches = zeros(max_read_len,2);
clear max_read_len
for i = 1:nreads
    read = deblank(reads(i,:));
    start = starts(i);
    len = length(read);
%     if (len >= 100 && start + len - 1 <= length(reference))
        subref = reference(start:start + len - 1);
        read = read(1:length(subref));
        len = length(read);
        these_mismatches = (read ~= subref & read ~= '-')';
        mismatches(1:len,1) = mismatches(1:len,1) + these_mismatches;
        mismatches(1:len,2) = mismatches(1:len,2) + ones(len,1);
%     end
end

clear nreads read start len subref these_mismatches

mean_mismatches = mismatches(:,1) ./ mismatches(:,2);
mean_mismatches = mean_mismatches(1:120);
h = figure(1);
plot(mean_mismatches,'k.');
ylim([0 1.1*max(mean_mismatches)]);
xlabel('position along read');
ylabel('fraction of mismatches to reference over all reads');
% print(h,'-deps','mismatches.eps')



% mean_mismatches = mismatches(:,1) ./ mismatches(:,2);
% mean_mismatches = mean_mismatches(1:120);
% h = figure(1);
% plot(mean_mismatches,'k.');
% ylim([0 1.1*max(mean_mismatches)]);
% xlabel('position along read');
% ylabel('fraction of mismatches to reference over all reads');


% h = figure(1);
% plot(mean_mismatches,'r.');
% hold on
% plot(mean_mismatches2,'b.');
% ylim([0 1.1*max(mean_mismatches)]);
% xlabel('position along read');
% ylabel('fraction of mismatches to reference over all reads');
% legend('cow pox', 'e. coli')
% legend('boxoff')