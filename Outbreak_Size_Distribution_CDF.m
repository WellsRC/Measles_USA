function [cdf_x] = Outbreak_Size_Distribution_CDF(x,k_nbin,p)

xv=[1:x];
log_pdf=gammaln(k_nbin+xv)-gammaln(k_nbin)-gammaln(xv+1)+k_nbin.*log(p)+xv.*log(1-p);

cdf_x=sum(exp(log_pdf));
end

