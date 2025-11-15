function [Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(F_NB,Max_Outbreak,p_c,mu_c,k_nbin,Reff,k_mealses,NS,Imported_Case)
rng(20251009)
r_z_samp=rand(size(p_c,1),NS);
r_os_samp=rand(size(p_c,1),NS);
Total_Cases_County=Imported_Case;
for ss=1:size(p_c,1)
    r_z=r_z_samp(ss,:);
    if(Reff(ss)>1)
        z_nb=F_NB(log10(k_nbin),log10(max(mu_c(ss),1)))';
        x0=1./(1+exp(-z_nb));
        % pd = makedist('NegativeBinomial','R',k_nbin,'P',p_nb);
        % pd = truncate(pd,1,inf);
        x=[1:Max_Outbreak(ss)];
        [lb,fval]=fmincon(@(z) (sum(x.*(logncdf(x,z,2)-logncdf(x-1,z,2)))./(sum(logncdf(x,z,2)-logncdf(x-1,z,2)))-mu_c(ss)).^2,0);
        p_nb=surrogateopt(@(z)  (sum(x.*(logncdf(x,z,2)-logncdf(x-1,z,2)))./(sum(logncdf(x,z,2)-logncdf(x-1,z,2)))-mu_c(ss)).^2,lb,lb+1);
        cdf_v=logncdf(x,p_nb,2);
        cdf_v=cdf_v./cdf_v(end);
        os=zeros(sum(r_z>p_c(ss,:)),1);
        temp_r=r_os_samp(ss,r_z>p_c(ss,:));
        for jj=1:length(os)
            fx=find(temp_r(jj)<=cdf_v,1,'first');
            os(jj)=x(fx);
        end
    else
        os=zeros(sum(r_z>p_c(ss,:)),1);
        temp_r=r_os_samp(ss,r_z>p_c(ss,:));
        for jj=1:length(os)
            cc=1;
            [pdf_0] = Chain_Size_Distribution(cc,Reff(ss),k_mealses);

            r=pdf_0+(1-pdf_0).*temp_r(jj);
            while (pdf_0<r && cc<min(101,Max_Outbreak(ss)))
                cc=cc+1;
                pdf_0 = pdf_0+Chain_Size_Distribution(cc,Reff(ss),k_mealses);
            end
            os(jj)=cc; % Need to remove the introductory seed for the chain as we classified it as an import in the likelihood
        end
    end
    Total_Cases_County(ss,r_z>p_c(ss,:))=Total_Cases_County(ss,r_z>p_c(ss,:))+os';
end
end