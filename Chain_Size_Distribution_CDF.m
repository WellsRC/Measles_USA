function [cdf_x] = Chain_Size_Distribution_CDF(x,R_eff,k)

xv=[1:x];
log_pdf=gammaln(k.*xv+xv-1)-gammaln(k.*xv)-gammaln(xv+1)+(xv-1).*log(R_eff./k)-(k.*xv+xv-1).*log(1+R_eff./k);

cdf_x=sum(exp(log_pdf));
end

