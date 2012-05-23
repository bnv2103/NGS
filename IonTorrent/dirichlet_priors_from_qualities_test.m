% get MLE for Dirichlet priors

addpath(genpath('lightspeed'));
addpath(genpath('fastfit'));

disp('genome wide allele fractions:')
c = [.9
     .05
     .03
     .02]

disp('depth:')
d = 1000

disp('pileup counts:')
n = [850
      80
      40
      30]

% specify true dirichlet parameters, and draw quality scores from that

disp('drawing quality scores from Dirichlet distributions...')
% concentration parameters for the four models
s_true = {10 10 10 10};

% means for the four models, random about genome-wide fractions
for i = 1:4
    m_true{i} = c + .1*randn(4, 1).*c;
    m_true{i} = m_true{i}./sum(m_true{i});
    assert(min(m_true{i}) >= 0);
end

% sample from these distributions to generate q data
for i = 1:4
    D_true{i} = dirichlet_sample(s_true{i}*m_true{i}, n(i));
end

% functions for phred scaling probability
phred = @ (p) -10*log10(p) + 33;
unphred = @ (q) 10.^(-(q-33)./10);

% quality scores
for i = 1:4
    q{i} = phred(1 - D_true{i}(i, :));
end

disp('generating data D using drawn quality scores and genome-wide fractions to fill in unknowns...')
      
D{1} = [1 - unphred(q{1})
        unphred(q{1})*c(2)/(1 - c(1))
        unphred(q{1})*c(3)/(1 - c(1))
        unphred(q{1})*c(4)/(1 - c(1))];
    
D{2} = [unphred(q{2})*c(1)/(1 - c(2))
        1 - unphred(q{2})
        unphred(q{2})*c(3)/(1 - c(2))
        unphred(q{2})*c(4)/(1 - c(2))];

D{3} = [unphred(q{3})*c(1)/(1 - c(3))
        unphred(q{3})*c(2)/(1 - c(3))
        1 - unphred(q{3})
        unphred(q{3})*c(4)/(1 - c(3))];

D{4} = [unphred(q{4})*c(1)/(1 - c(4))
        unphred(q{4})*c(2)/(1 - c(4))
        unphred(q{4})*c(3)/(1 - c(4))
        1 - unphred(q{4})];

% each of these can be used for MLE of the respective prior

disp('recovering Dirichlet parameters from the D by Maximum likelihood...')

% set this to 'true' to test with D_true instead of D, which is made with genome-wide ratios
SANITY_CHECK = false;

for i = 1:4
    alpha_guess = 0.5*c';
    if SANITY_CHECK
        alpha = dirichlet_fit(D_true{i}', alpha_guess);
    else
        alpha = dirichlet_fit(D{i}', alpha_guess);
    end
    s{i} = sum(alpha);
    m{i} = alpha'/s{i};
end

disp(' ')
disp('results:')
disp(' ')
disp('true Dirichlet parameters:')
disp('s = ')
disp(cell2mat(s_true))
disp('m = ')
disp(cell2mat(m_true))
if SANITY_CHECK
    disp('computed Dirichlet parameters from synthetic D, NOT using genome-wide ratios:')
else
    disp('computed Dirichlet parameters:')
end
disp('s = ')
disp(cell2mat(s))
disp('m = ')
disp(cell2mat(m))
