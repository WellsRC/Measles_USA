function [Outbreak_County]=Monte_Carlo_Outbreak_County_Fitting(F_NB,Max_Outbreak,p_c,mu_c,k_nbin,Reff,k_mealses,r_samp_pc,r_samp_outbreak)
Outbreak_County=zeros(size(r_samp_pc));
for ss=1:length(p_c)
    r_z=r_samp_pc(ss,:);
    r_outbreak=r_samp_outbreak(ss,:);
    r_outbreak=r_outbreak(r_z>p_c(ss));
    if(Reff(ss)>1)
        z_nb=F_NB(log10(k_nbin),log10(max(mu_c(ss),1)))';
        p_nb=1./(1+exp(-z_nb));
        pd = makedist('NegativeBinomial','R',k_nbin,'P',p_nb);
        pd = truncate(pd,1,inf);
        x=[1:Max_Outbreak(ss)];
        cdf_v=cdf(pd,x);
        cdf_v=cdf_v./cdf_v(end);
        os=zeros(sum(r_z>p_c(ss)),1);
        for jj=1:length(os)
            fx=find(r_outbreak(jj)<=cdf_v,1,'first');
            os(jj)=x(fx);
        end
    else
        os=zeros(sum(r_z>p_c(ss)),1);
        for jj=1:length(os)
            cc=1;
            [pdf_0] = Chain_Size_Distribution(cc,Reff(ss),k_mealses);

            r=pdf_0+(1-pdf_0).*r_outbreak(jj);
            while (pdf_0<r && cc<min(101,Max_Outbreak(ss)))
                cc=cc+1;
                pdf_0 = pdf_0+Chain_Size_Distribution(cc,Reff(ss),k_mealses);
            end
            os(jj)=cc; % Need to remove the introductory seed for the chain as we classified it as an import in the likelihood
        end
    end
    Outbreak_County(ss,r_z>p_c(ss))=os;
end
end