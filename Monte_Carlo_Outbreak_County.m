function [Outbreak_County]=Monte_Carlo_Outbreak_County(p_c,mu_c,k_nbin,Reff,k_mealses,NS)
rng(20251009)
Outbreak_County=zeros(length(p_c),NS);
for ss=1:length(p_c)
    r_z=rand(NS,1);
    if(Reff(ss)>1)
        os=nbinrnd(k_nbin,k_nbin./(k_nbin+mu_c(ss)),sum(r_z>p_c(ss)),1);
    else
        os=zeros(sum(r_z>p_c(ss)),1);
        for jj=1:length(os)
            cc=1;
            r=rand(1);
            [pdf_0] = Chain_Size_Distribution(cc,Reff(ss),k_mealses);
            while (pdf_0<r & cc<100)
                cc=cc+1;
                pdf_0 = pdf_0+Chain_Size_Distribution(cc,Reff(ss),k_mealses);
            end
            os(jj)=cc-1; % Need to remove the introductory seed for the chain as we classified it as an import in the likelihood
        end
    end
    Outbreak_County(ss,r_z>p_c(ss))=os;
end
end