[FileName,PathName,FilterIndex] = uigetfile('*mpileup.cts', 'MultiSelect', 'on');

for i = 1:numel(FileName)
    thisfile = FileName{i};
    data = importdata([PathName,thisfile],'\t',1);
    data = data.data;

    % bar(data(:,1)-178936093, data(:,5:7),'stacked')
    % hold on
    % plot(data(:,1)-178936093, .01*data(:,2))
    % xlabel('base distance from 178936093')
    % ylabel('depth of coverage')
    % legend('allele 1', 'allele 2', 'allele 3', '1% coverage')
    % print('-f1','-depsc','parsepileup.eps')   

    ind = find(data(:,2) >= 0);

    figure(1)
    bar(data(ind,1), data(ind,5:7),'stacked')
    hold on
    plot(data(ind,1), .01*data(ind,2))
    xlabel('base position')
    ylabel('depth of coverage')
    % ylim([0 50])
    legend('allele 1', 'allele 2', 'allele 3', '1% coverage')
    slashes = find(PathName == '/');
    gene = PathName(slashes(end-1)+1:end-1);
    title(['gene: ', gene, ', ', 'amplicon: ', thisfile(1:end-12)])
    hold off
    print('-f1','-dpdf',[PathName,thisfile,'.pdf'])
    close 1
end