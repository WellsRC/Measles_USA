function [x_values,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_baseline,pd_reduction_1,pd_reduction_2,pd_reduction_3)

    min_x=floor(min([icdf(pd_baseline,0.001) icdf(pd_reduction_1,0.001) icdf(pd_reduction_2,0.001) icdf(pd_reduction_3,0.001)]));
    max_x=ceil(max([icdf(pd_baseline,0.999) icdf(pd_reduction_1,0.9) icdf(pd_reduction_2,0.9) icdf(pd_reduction_3,0.9)]));
    
    dx=[1 5 10 25 50 100 250 500 1000 5000 10^4 2.5.*10^4 5.*10^4 10^5 2.5.*10^5 5.*10^5 10^6 2.5.*10^6 5.*10^6 10^7 2.5.*10^7 5.*10^7 10^8 2.5.*10^8 5.*10^8 10^9 2.5.*10^9 5.*10^9];
    Num_indx=zeros(length(dx),1);
    for jj=1:length(dx)
        Num_indx(jj)=length((min_x-rem(min_x,dx(jj))):dx(jj):(ceil(max_x./dx(jj)).*dx(jj)));
    end
    fval=(Num_indx-11).^2;
    dx=max(dx(fval==min(fval)));
    if(min_x<dx)
        min_x=0;
    else
        min_x=min_x-rem(min_x,dx);
    end

    max_x=ceil(max_x./dx).*dx;

    x_values=linspace(min_x,max_x,10^4);
end