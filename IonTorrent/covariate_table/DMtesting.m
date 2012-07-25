path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/lightspeed', path);
path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/fastfit', path);

alphai = 1000*[.994 .002 .002 .002];

dat1 = load('cellcounts/A_0');
alpha1 = polya_fit(dat1, alphai);
dat1 = dat1(find(sum(dat1, 2) > 1000), :);

dat2 = load('cellcounts/A_4_0');
alpha2 = polya_fit(dat2, alphai);
dat2 = dat2(find(sum(dat2, 2) > 1000), :);



% mismatch prob dist

dat1_collapsed = sum(dat1(:, 2:4), 2)./sum(dat1, 2);
dat2_collapsed = sum(dat2(:, 2:4), 2)./sum(dat2, 2);

binsize = .0005;
binedges = 0:binsize:.05;
cts1 = histc(dat1_collapsed, binedges)/(binsize*numel(dat1_collapsed));
cts2 = histc(dat2_collapsed, binedges)/(binsize*numel(dat2_collapsed));

bar(binedges, [cts1 cts2])
hold on
fit1 = betapdf(binedges, sum(alpha1(2:4)), alpha1(1));
plot(binedges, fit1, 'b', 'LineWidth',2)
fit2 = betapdf(binedges, sum(alpha2(2:4)), alpha2(1));
plot(binedges, fit2, 'r', 'LineWidth',2)
hold off
xlim([0 .05])
legend('A-0', 'A-4-0', 'fit1', 'fit2')