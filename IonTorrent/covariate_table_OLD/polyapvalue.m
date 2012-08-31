% get polya p value
function pval = polyapvalue(cts, p, ref_cat)

% depth
d = sum(cts);

% permute indices to put reference category first
cts = [cts(ref_cat) cts(1:ref_cat - 1) cts(ref_cat + 1:end)];
p = [p(:, ref_cat) p(:, 1:ref_cat - 1) p(:, ref_cat + 1:end)];

alpha = dirichlet_fit_newton(p);

% easier to compute complement (less iteration)
pval_comp = 0;
for n1 = (cts(1) + 1):d
    for n2 = 0:(d - n1)
        for n3 = 0:(d - n2 - n1)
            n = [n1 n2 n3 (d - n3 - n2 - n1)];
            pval_comp = pval_comp + sum(exp(polya_logProb(alpha, n)));
        end
    end
end

pval = 1 - pval_comp;

end