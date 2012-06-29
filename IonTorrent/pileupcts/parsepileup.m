sample_folder = uigetdir;

gene_folders = dir(sample_folder);

for g = 1:numel(gene_folders)

if ~gene_folders(g).isdir || gene_folders(g).name(1) == '.'; continue; end

cts_files = dir([sample_folder, '/', gene_folders(g).name, '/*cts']);

% [FileName,PathName,FilterIndex] = uigetfile('*mpileup.cts', 'MultiSelect', 'on');

namps = numel(cts_files);

%namps = numel(FileName);

figure(1);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 30 15]);

plt = 1;

for i = 1:namps
    %thisfile = FileName{i};
    thisfile = cts_files(i).name;
    PathName = [sample_folder, '/', gene_folders(g).name, '/'];
    %data = importdata([PathName,thisfile],'\t',1);
    data = importdata([PathName, thisfile],'\t',1);
    if ~isstruct(data); continue; end
    data = data.data;

    % we also need to grab the chromosome name from the mpileup file
    %mpfile = [PathName, thisfile(1:end-4)];
    mpfile = [PathName, thisfile(1:end-4)];
    fid = fopen(mpfile, 'r');
    chr_str = strtok(fgetl(fid));
    fclose(fid);

    % threshold on coverage
    ind = find(data(:,2) >= 0);
    % candidate mutation sites have top nonreference exceed 1% coverage
    candidate = find(data(ind,5) > .01*data(ind,2));

    subplot(1, namps, plt) 
    plt = plt + 1;
    bar(data(ind,1), data(ind,5:7),'stacked')
    for c=candidate
        text(data(ind(c),1), 1.01*sum(data(ind(c),5:7), 2), num2str(data(ind(c),1)), 'FontSize', 8, 'HorizontalAlignment', 'center');
    end

    hold on
    plot(data(ind,1), .01*data(ind,2))
    hold off
    if i == 1
        xlabel(['base position in ', chr_str])
        ylabel('depth of coverage')
    end
    set(gca,'YScale','log')
    set(gca,'XTick',data(ind([1 end]), 1))
    set(gca,'XTickLabel',{num2str(data(ind(1), 1)), num2str(data(ind(end), 1))})
    xlim([data(ind([1])) data(ind([end]))])
    %ylim([0 1.1*max(sum(data(ind,5:7), 2))])
    ylim([0 10000])
    pbaspect([1 2 1])
    %axis square
    % ylim([0 50])
    %legend('alt. allele 1', 'alt. allele 2', 'alt. allele 3', '1% total coverage')
    %legend('boxoff')
    slashes = find(PathName == '/');
    sample = strtok(PathName(slashes(end-2)+1:slashes(end-1)-1), '.');
    gene = PathName(slashes(end-1)+1:end-1);
    title(thisfile(1:end-12))
end

%ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%text(0.5, 1, ['sample: ', sample, ', ', 'gene: ', gene], 'HorizontalAlignment' ,'center','VerticalAlignment', 'top')

%print('-f1','-depsc', [PathName, sample, '_', gene, '.eps'])
print('-f1','-depsc', [PathName, sample, '_', gene, '.eps'])
%    close all 
%end

%close all


end
