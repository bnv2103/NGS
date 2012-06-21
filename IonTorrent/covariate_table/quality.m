% plot figures from output of quality.py

base = input('input base name: ', 's');

disp('loading read length data')
L = importdata([basename '.lengths']);
mu_L = mean(L);
sigma_L = std(L);

disp('loading forward read data...')
disp(' A...')
Af = importdata([basename '.forwardA.cor']);
Af = Af(find(Af(:,3) >= 0 & Af(:,4) - mu_L <= sigma_L),:);
Af(:,1) = Af(:,1) - 33;
disp(' T...')
Tf = importdata([basename '.forwardT.cor']);
Tf = Tf(find(Tf(:,3) >= 0 & Tf(:,4) - mu_L <= sigma_L),:);
Tf(:,1) = Tf(:,1) - 33;
disp(' C...')
Cf = importdata([basename '.forwardC.cor']);
Cf = Cf(find(Cf(:,3) >= 0 & Cf(:,4) - mu_L <= sigma_L),:);
Cf(:,1) = Cf(:,1) - 33;
disp(' G...')
Gf = importdata([basename '.forwardG.cor']);
Gf = Gf(find(Gf(:,3) >= 0 & Gf(:,4) - mu_L <= sigma_L),:);
Gf(:,1) = Gf(:,1) - 33;

disp('loading reverse read data...')
disp(' A...')
Ar = importdata([basename '.reverseA.cor']);
Ar = Ar(find(Ar(:,3) >= 0 & Ar(:,4) - mu_L <= sigma_L),:);
Ar(:,1) = Ar(:,1) - 33;
disp(' T...')
Tr = importdata([basename '.reverseT.cor']);
Tr = Tr(find(Tr(:,3) >= 0 & Tr(:,4) - mu_L <= sigma_L),:);
Tr(:,1) = Tr(:,1) - 33;
disp(' C...')
Cr = importdata([basename '.reverseC.cor']);
Cr = Cr(find(Cr(:,3) >= 0 & Cr(:,4) - mu_L <= sigma_L),:);
Cr(:,1) = Cr(:,1) - 33;
disp(' G...')
Gr = importdata([basename '.reverseG.cor']);
Gr = Gr(find(Gr(:,3) >= 0 & Gr(:,4) - mu_L <= sigma_L),:);
Gr(:,1) = Gr(:,1) - 33;

disp('plotting...')

set(0,'DefaultFigureWindowStyle','docked')


figure(1)
bins = 0:max(L);
vals = histc(L, bins);
vals = vals/sum(vals);
pass = bins >= mu_L - sigma_L & bins <= mu_L + sigma_L;
bar(bins(~pass), vals(~pass), 'r')
hold on
bar(bins(pass), vals(pass), 'g')
hold off
axis([0 max(L) 0 1.1*max(vals)])
xlabel('read length')
legend('filtered', 'passed filter')

figure(2)
bins = 4:29;
N = [histc(Af(:,1), bins) histc(Tf(:,1), bins), histc(Cf(:,1), bins) histc(Gf(:,1), bins)];
N = N./repmat(sum(N,1),size(N,1),1);
semilogy(10.^(-bins/10), N, '->', 'LineWidth', 2, 'MarkerSize', 10)
hold on
N = [histc(Ar(:,1), bins) histc(Tr(:,1), bins), histc(Cr(:,1), bins) histc(Gr(:,1), bins)];
N = N./repmat(sum(N,1),size(N,1),1);
semilogy(10.^(-bins/10), N, '-<', 'LineWidth', 2, 'MarkerSize', 10)
hold off
xlabel('p_{mismatch} = 10^{ -q / 10}')
legend('forward A', 'forward T', 'forward C', 'forward G', ...
       'reverse A', 'reverse T', 'reverse C', 'reverse G')
axis([0 .33 0 1.1*max(N(:))])


levels = log(logspace(log(0.001),0));

figure(3)
q = 0:35;
d = 0:20;
N = hist3([Af(:,1) Af(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'b')
hold on
N = hist3([Tf(:,1) Tf(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'g')
N = hist3([Cf(:,1) Cf(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'r')
N = hist3([Gf(:,1) Gf(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'c')
hold off
title('forward reads')
xlabel('nucleotide distance from nearest homopolymer 3'' boundary')
ylabel('quality score')
legend('A', 'T', 'C', 'G')

figure(4)
N = hist3([Ar(:,1) Ar(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'b')
hold on
N = hist3([Tr(:,1) Tr(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'g')
N = hist3([Cr(:,1) Cr(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'r')
N = hist3([Gr(:,1) Gr(:,3)], 'Edges', {q, d});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(d, q, N, levels, 'c')
hold off
title('reverse reads')
xlabel('nucleotide distance from nearest homopolymer 3'' boundary')
ylabel('quality score')
legend('A', 'T', 'C', 'G')

figure(5)
p = 0:mu_L + sigma_L;
N = hist3([Af(:,1) Af(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'b')
hold on
N = hist3([Tf(:,1) Tf(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'g')
N = hist3([Cf(:,1) Cf(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'r')
N = hist3([Gf(:,1) Gf(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'c')
hold off 
title('forward reads')
xlabel('nucleotide position from 5'' end')
ylabel('quality score')
legend('A', 'T', 'C', 'G')

figure(6)
N = hist3([Ar(:,1) Ar(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'b')
hold on
N = hist3([Tr(:,1) Tr(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'g')
N = hist3([Cr(:,1) Cr(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'r')
N = hist3([Gr(:,1) Gr(:,2)], 'Edges', {q, p});
N = log(N./repmat(sum(N,1),size(N,1),1));
contour(p, q, N, levels, 'c')
hold off
title('reverse reads')
xlabel('nucleotide position from 5'' end')
ylabel('quality score')
legend('A', 'T', 'C', 'G')

disp('saving plots...')
print('-f1','-dpdf',[basename '.plot1.pdf'])
close(1)
print('-f2','-dpdf',[basename '.plot2.pdf'])
close(2)
print('-f3','-dpdf',[basename '.plot3.pdf'])
close(3)
print('-f4','-dpdf',[basename '.plot4.pdf'])
close(4)
print('-f5','-dpdf',[basename '.plot5.pdf'])
close(5)
print('-f6','-dpdf',[basename '.plot5.pdf'])
close(6)
