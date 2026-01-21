function [mle,lb_hdi,ub_hdi]=Estimate_HDI(p,hdi)
    c=linspace(0,icdf(p,0.95),1001);
    x0=c(pdf(p,c)==max(pdf(p,c)));
    mle=fmincon(@(z)-log(pdf(p,z)),x0);

    lb0=(mle-icdf(p,0))./5;
    bnd=fmincon(@(z) (log(pdf(p,mle-z(1)))-log(pdf(p,mle+z(2)))).^2+(100.*cdf(p,mle+z(2))-100.*cdf(p,mle-z(1))-100.*hdi).^2,[lb0 mle./5],[],[],[],[],[0 0],[mle inf]);

    lb_hdi=mle-bnd(1);
    ub_hdi=mle+bnd(2);
end