function [pdf_x] = Chain_Size_Distribution(x,R_eff,k)


log_pdf=gammaln(k.*x+x-1)-gammaln(k.*x)-gammaln(x+1)+(x-1).*log(R_eff./k)-(k.*x+x-1).*log(1+R_eff./k);

pdf_x=exp(log_pdf);
end

