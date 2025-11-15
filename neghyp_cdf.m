function Y = neghyp_cdf(x,N,K,R)

xv=[0:x];
Y=min(sum(neghyp_pdf(xv,N,K,R)),1);

end