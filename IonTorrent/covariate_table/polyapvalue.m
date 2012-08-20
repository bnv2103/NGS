% get polya p value given counts and alpha
function [pval, alt] = polyapvalue(cts, alpha, ref)

% identify top nonref allele. Note, previous germline filtering step garantees that fraction of reference is >= 60%
[cts, sortidx] = sort(cts, 'descend');
% must be a dominant nonref
if ~(cts(2) > cts(3))
    pval = NaN;
    alt = '';
    return
end

alpha = alpha(sortidx);
alleles = 'ACGT';
alt = alleles(sortidx(2));

% depth
d = sum(cts);

pval = 0;
for n1 = cts(1):-1:0
    n2 = cts(2) + (cts(1) - n1);
    for n3 = 0:(d - n1 - n2)
        n4 = d - n3 - n2 - n1;
        n = [n1 n2 n3 n4];

%       note: Minka's Polya is for the emmision case (order specific)
%       we have to multiply by a multinomial coefficient
%       overflows for large counts, Inf/Inf = NaN, so
%       we're just ignoring those terms, could be problematic.

%       we can use Stirling's approximation for the multinomial coeff, we
%       know that factorial overflows for arguments over 170
        if d > 170
            log_d_fact = d*log(d) - d;
        else
            log_d_fact = log(factorial(d));
        end
        log_n_fact = zeros(4, 1);
        for ni = 1:4
            if n(ni) > 170
                log_n_fact(ni) = n(ni)*log(n(ni)) - n(ni);
            else
                log_n_fact(ni) = log(factorial(n(ni)));
            end
        end
        log_multicoef = log_d_fact - sum(log_n_fact);

        pval = pval + exp(log_multicoef + polya_logProb(alpha, n));   
    end
end

end
