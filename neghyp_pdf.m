function y = neghyp_pdf(x,N,K,R)

try
    temp_1=gammaln(x+R)-gammaln(x+1)-gammaln(R);
    temp_2=gammaln(N-R-x+1)-gammaln(K-x+1)-gammaln(N-R-K+1);
    temp_3=gammaln(N+1)-gammaln(K+1)-gammaln(N-K+1);
    
    y=exp(temp_1+temp_2-temp_3);
catch
    y=NaN;
end

end