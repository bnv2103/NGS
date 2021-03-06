% get MLE for Dirichlet prior, null model

addpath(genpath('lightspeed'));
addpath(genpath('fastfit'));

[file, path] = uigetfile('*.puq');

fin = fopen([path '/' file]);

puq = textscan(fin, '%c %s %s %s %s %s %s %s %s', ...
                    'Delimiter', '\t', 'BufSize', 10000);

nsites = numel(puq{1});

% genome-wide counts
N = zeros(4, 4, 2); 

for site = 1:nsites
    for forward = [true false]
        if forward
            dir = 0;
            row = 1;
        else
            dir = 4; % index offset to get reverse read data
            row = 2;
        end

        reference = puq{1}(site);
        refind = find(['A' 'T' 'C' 'G'] == reference);

        qA = cell2mat(puq{2 + dir}(site));
        qT = cell2mat(puq{3 + dir}(site));
        qC = cell2mat(puq{4 + dir}(site));
        qG = cell2mat(puq{5 + dir}(site));

        % counts
        N(refind, :, row) = N(refind, :, row) + [length(qA) length(qT) length(qC) length(qG)];
    end
end

% eventually we will loop over all pileups, forward and reverse
for site = 1:nsites
    for forward = [true false]
        if forward
            dir = 0;
            row = 1;
        else
            dir = 4; % index offset to get reverse read data
            row = 2;
        end

        
        
        reference = puq{1}(site);
        refind = find(['A' 'T' 'C' 'G'] == reference);

        qA = cell2mat(puq{2 + dir}(site));
        qT = cell2mat(puq{3 + dir}(site));
        qC = cell2mat(puq{4 + dir}(site));
        qG = cell2mat(puq{5 + dir}(site));

        % counts
        n = [length(qA) length(qT) length(qC) length(qG)];
        
        if sum(n) < 1000
            continue
        end
        

        
        [n ni] = sort(n, 'descend');
        
        
        % kludge
        if find(['A' 'T' 'C' 'G'] == reference) ~= ni(1)
            continue
        end        
        
        
        assert(find(['A' 'T' 'C' 'G'] == reference) == ni(1))

        % genome wide allele fractions for this reference base
        Ntot = sum(N(refind, :, row));
        for i = 1:4
            c(i) = N(refind, ni(i), row)/Ntot;
        end

        % depth
        d = sum(n);

        q{1} = double(cell2mat(puq{1 + ni(1) + dir}(site)));
        q{2} = double(cell2mat(puq{1 + ni(2) + dir}(site)));
        q{3} = double(cell2mat(puq{1 + ni(3) + dir}(site)));
        q{4} = double(cell2mat(puq{1 + ni(4) + dir}(site)));

        % functions for phred scaling probability
        phred = @ (p) -10*log10(p) + 33;
        unphred = @ (q) 10.^(-(q-33)./10);

%         disp('generating data D using quality scores and genome-wide fractions to fill in unknowns...')

        % null model m1
        % reads called as a1
        D(:, 1:numel(q{1})) = [1 - unphred(q{1})
                               unphred(q{1})*c(2)/(1 - c(1))
                               unphred(q{1})*c(3)/(1 - c(1))
                               unphred(q{1})*c(4)/(1 - c(1))];

        % reads called as a2
        if n(2)
            D(:, end + 1:end + numel(q{2})) = [(1 - unphred(q{2}))*c(1)/(1 - c(2))
                                               unphred(q{2})
                                               (1 - unphred(q{2}))*c(3)/(1 - c(2))
                                               (1 - unphred(q{2}))*c(4)/(1 - c(2))]; 
        end

        % reads called as a3
        if n(3)
            D(:, end + 1:end + numel(q{3})) = [(1 - unphred(q{3}))*c(1)/(1 - c(3))
                                               (1 - unphred(q{3}))*c(2)/(1 - c(3))
                                               unphred(q{3})
                                               (1 - unphred(q{3}))*c(4)/(1 - c(3))]; 
        end

        % reads called as a4
        if n(4)
            D(:, end + 1:end + numel(q{4})) = [(1 - unphred(q{4}))*c(1)/(1 - c(4))
                                                  (1 - unphred(q{4}))*c(2)/(1 - c(4))
                                                  (1 - unphred(q{4}))*c(3)/(1 - c(4))
                                                  unphred(q{4})];
        end
        % MLE

%         disp('recovering Dirichlet parameters from D by Maximum likelihood...')

        alpha_guess = 0.5*c;
        alpha = dirichlet_fit(D', alpha_guess);
%         s = sum(alpha);
%         m = alpha'/s;

%         disp(' ')
%         disp('results:')
%         disp(' ')
%         disp('MLE computed Dirichlet parameters:')
%         disp('s = ')
%         disp(s)
%         disp('m = ')
%         disp(m)

        % matrix of all possible count vectors consistent with n(1) and d
        n_prime = [];
        for n_2 = 0:d - n(1)
            for n_3 = 0:d - n(1) - n_2 - 1
                n_prime(end + 1, :) = [n(1) n_2 n_3 (d - n(1) - n_2 - n_3)];
            end
        end

        p_value(site) = sum(exp(polya_logProb(alpha, n_prime)));

%         p_value(site) = exp(polya_logProb(alpha, n));



        % fitting marginal betas and getting alpha from that 
%         shape{1} = betafit(1-unphred(q{1}));
%         shape{2} = betafit(unphred(q{2}));
%         shape{3} = betafit(unphred(q{3}));
%         shape{4} = betafit(unphred(q{4}));
% 
%         b = [shape{1}(1)
%              shape{2}(1)
%              shape{3}(1)
%              shape{4}(1)
%              shape{1}(2)
%              shape{2}(2)
%              shape{3}(2)
%              shape{4}(2)];
% 
%         A = [eye(4); ~eye(4)];
% 
%         % alpha = lsqlin(A, b, [], [], [], [], zeros(4, 1), [], s*m, optimset('MaxIter', Inf, 'TolFun', 1.0e-12, 'Display', 'on'));
%         problem.C = A;
%         problem.d = b;
%         problem.x0 = s*m;
%         problem.solver = 'lsqnonneg';
%         problem.options = optimset('Display', 'final', 'TolX', 1.0e-9);
% 
%         alpha = lsqnonneg(problem);
% 
%         disp('Dirichlet parameters computed from marginal beta MLEs:')
%         s_dumb = sum(alpha);
%         m_dumb = alpha/s_dumb;
% 
%         disp('s = ')
%         disp(s_dumb)
%         disp('m = ')
%         disp(m_dumb)        




%     pause



    end
end

fclose(fin);

plot(p_value)
