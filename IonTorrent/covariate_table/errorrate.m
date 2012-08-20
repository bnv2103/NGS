% collapse Dirichlets and plot error rate

clear all; close all; clc

x = 0:.0001:1;

nucs = 'ACGT';

for b = 1:4
    subplot(2, 2, b)
    Dvector = load([nucs(b) '_0.alpha']);
    beta = Dvector(b);
    alpha = sum(Dvector) - Dvector(b);
    plot(x, betapdf(x, alpha, beta)/length(x), 'k')
    hold on
                
    files = dir([nucs(b) '_*_0.alpha']);
    nfiles = numel(files);
    for f = 1:nfiles
        Dvector = load(files(f).name);
        beta = Dvector(b);
        alpha = sum(Dvector) - Dvector(b);
        dat(:, f) = betapdf(x, alpha, beta)'/length(x);
        files(f).legend = files(f).name(3:end-6);
    end
    plot(x, dat)
    
    if b == 1
        legend('< 3', '3', '> 3')
    end
    
    files = dir([nucs(b) '*_*_1.alpha']);
    nfiles = numel(files);
    for f = 1:nfiles
        Dvector = load(files(f).name);
        beta = Dvector(b);
        alpha = sum(Dvector) - Dvector(b);
        dat(:, f) = betapdf(x, alpha, beta)'/length(x);
        files(f).legend = files(f).name(3:end-6);
    end
    plot(x, dat, '--')
    xlim([0 .02])
    xlabel('error probability')
    title(nucs(b))
    hold off
end

% legend(files.legend)

pause
legend boxoff
print(1, '-depsc', 'errorrate.eps')