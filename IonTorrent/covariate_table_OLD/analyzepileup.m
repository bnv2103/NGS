% compute alpha vectors for pileup

addpath(genpath('../matlab_code/lightspeed'));
addpath(genpath('../matlab_code/fastfit'));

matdir = uigetdir;
matfiles = dir([matdir, '/*.mat']);

catorder = ['A' 'C' 'G' 'T'];

% thresh = 1000;

pos_axis = [];
pval_for = [];
pval_rev = [];

for i = 1:numel(matfiles)
    
    load([matdir, '/', matfiles(i).name]);
    C = textscan(matfiles(i).name, '%[^_]_%u_%c%c.mat');
    chr = C{1};
    pos = C{2};
    ref_for = C{3};
    ref_rev = C{4};
    ref_for_cat = find(catorder == ref_for);
    ref_rev_cat = find(catorder == ref_rev);
    
    totaltablecounts_for = sum(tablecounts_for, 2);
    totaltablecounts_rev = sum(tablecounts_rev, 2);
    
%     % only do something if for and rev covariate counts are high enough
%     if numel(totaltablecounts_for) < thresh || ...
%        numel(totaltablecounts_rev) < thresh
%         continue
%     end

%     % average phat and standard multinomial
%     phat_for_mean = sum(double(tablecounts_for));
%     phat_for_mean = phat_for_mean/sum(phat_for_mean);
%     phat_rev_mean = sum(double(tablecounts_rev));
%     phat_rev_mean = phat_rev_mean/sum(phat_rev_mean);      
%     pval_for(end + 1) = mnpvalue(double(counts_for'), phat_for_mean, ref_for_cat);
%     pval_rev(end + 1) = mnpvalue(double(counts_rev'), phat_rev_mean, ref_rev_cat);
    
    % polya
    phat_for = double(tablecounts_for)./repmat(totaltablecounts_for, 1, 4);
    phat_rev = double(tablecounts_rev)./repmat(totaltablecounts_rev, 1, 4);
    pval_for(end + 1) = polyapvalue(double(counts_for'), phat_for, ref_for_cat);
    pval_rev(end + 1) = polyapvalue(double(counts_rev'), phat_rev, ref_rev_cat);
    
    pos_axis(end + 1) = pos;
    
    clc
    disp('progress:')
    disp(i/numel(matfiles))
    
end

bar(pos_axis, [log10(pval_for); -log10(pval_rev)]')