function [pdf_x] = Chain_Size_Distribution(x,R_eff,k)

if(x<=100) % Truncated it at 100
    [cdf_x] = Chain_Size_Distribution_CDF(100,R_eff,k);
    log_pdf=gammaln(k.*x+x-1)-gammaln(k.*x)-gammaln(x+1)+(x-1).*log(R_eff./k)-(k.*x+x-1).*log(1+R_eff./k);
    
    pdf_x=exp(log_pdf)./cdf_x;
else
    pdf_x=0;
end
end

