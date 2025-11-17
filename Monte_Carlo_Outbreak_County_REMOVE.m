function [Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(Max_Outbreak,p_c,N_NHG,K_NHG,R_NHG,Reff,k_mealses,NS,Imported_Case)
rng(20251009)
r_z_samp=rand(size(p_c,1),NS);
r_os_samp=rand(size(p_c,1),NS);
Total_Cases_County=Imported_Case;
for ss=1:size(p_c,1)
    r_z=r_z_samp(ss,:);
    if(Reff(ss)>1)
        x=[1:Max_Outbreak(ss)]-1;
        pdf_v=neghyp_pdf(x,N_NHG(ss),K_NHG(ss),R_NHG);
        cdf_v=cumsum(pdf_v)./sum(pdf_v);
        os=zeros(sum(r_z>p_c(ss,:)),1);
        temp_r=r_os_samp(ss,r_z>p_c(ss,:));
        for jj=1:length(os)
            fx=find(temp_r(jj)<=cdf_v,1,'first');
            os(jj)=x(fx)+1; % ADD ONE AS we are using the negative hyper-geometric a the truncated 
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