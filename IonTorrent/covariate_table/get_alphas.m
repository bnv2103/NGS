function get_alphas(ctsfile, outfile)
% add Minka's code to the path
path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/lightspeed', path);
path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/fastfit', path);
% get the count files
dat = load(ctsfile);
% initial guess for alpha
alphai = 1000*[.002 .002 .002 .002];
[~, fname, ~] = fileparts(ctsfile);
alphai(find(fname(1) == 'ACGT')) = 0.994;
alpha = polya_fit(dat, alphai);
dlmwrite(outfile, alpha, '\t');
end
