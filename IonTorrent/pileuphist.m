% vuiimport('pileuphist*.out')v
% pause
% bar(data(:,1)-178936093, data(:,5:7),'stacked')
% hold on
% plot(data(:,1)-178936093, .01*data(:,2))
% xlabel('base distance from 178936093')
% ylabel('depth of coverage')
% legend('allele 1', 'allele 2', 'allele 3', '1% coverage')
% % print('-f1','-depsc','parsepileup.eps')

sample = importdata('parsepileup_negativecontrol.out','\t',3);
negativecontrol = importdata('parsepileup_sample.out','\t',3);

sample_ind = find(sample.data(:,2) >= 1);
negativecontrol_ind = find(negativecontrol.data(:,2) >= 1);

mismatchratio_sample = sample.data(sample_ind,3)./sample.data(sample_ind,2);
mismatchratio_negativecontrol = negativecontrol.data(negativecontrol_ind,3)./negativecontrol.data(negativecontrol_ind,2);

[n_s, xout] = hist(mismatchratio_sample, 50);
% [n_s, xout_s] = hist(sample.data(:,3)./sample.data(:,2), 50)
% [n_c, xout_c] = hist(negativecontrol.data(:,3)./negativecontrol.data(:,2), xout_s)

figure(1)
hist(mismatchratio_sample,xout)
hold on
hist(mismatchratio_negativecontrol,xout)
xlabel('mismatch fraction')
hold off
legend('sample','control')
xlim([0 1])
ylim([0 2000])
title('distribution of mismatch as a fraction of coverage depth')
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','r','EdgeColor','k','facealpha',0.5)
set(h(2),'FaceColor','b','EdgeColor','k','facealpha',0.8)


% print('-f1','-depsc','parsepileup.eps')


% figure(1)
% hist(mismatchratio_sample, xout)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
% hold on
% hist(mismatchratio_negativecontrol, xout)
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w','facealpha',0.75)
% xlabel('mismatch fraction')
% hold off
% legend('sample','control')