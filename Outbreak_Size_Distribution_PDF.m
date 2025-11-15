function [pdf_x] = Outbreak_Size_Distribution_PDF(x,k_nbin,p)

log_pdf=gammaln(k_nbin+x)-gammaln(k_nbin)-gammaln(x+1)+k_nbin.*log(p)+x.*log(1-p);

pdf_x=exp(log_pdf);
end

