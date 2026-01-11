function [Outbreak_County]=Monte_Carlo_Outbreak_County_Fitting(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,r_samp_pc,r_samp_outbreak)
Outbreak_County=zeros(size(r_samp_pc));
for ss=1:size(p_zero,1)
    d_zero=r_samp_pc(ss,:)-p_zero(ss,:);
    r_outbreak=r_samp_outbreak(ss,:);
    r_outbreak=r_outbreak(d_zero>0);
    if(Reff(ss)>1)
        x=[1:Max_Outbreak(ss)]-1;
        % x=[0:K_NHG(ss)];
        pdf_v=neghyp_pdf(x,N_NHG(ss),K_NHG(ss),R_NHG);
        cdf_v=cumsum(pdf_v)./sum(pdf_v);
        os=zeros(sum(d_zero>0),1);
        for jj=1:length(os)
            fx=find(r_outbreak(jj)<=cdf_v,1,'first');
            os(jj)=x(fx)+1; % ADD ONE AS we are using the negative hyper-geometric a the truncated 
        end
    else
        os=zeros(sum(d_zero>0),1);
        for jj=1:length(os)
            cc=1;
            [pdf_0] = Chain_Size_Distribution(cc,Reff(ss),k_mealses);

            r=pdf_0+(1-pdf_0).*r_outbreak(jj);
            while (pdf_0<r && cc<min(100,Max_Outbreak(ss)))
                cc=cc+1;
                pdf_0 = pdf_0+Chain_Size_Distribution(cc,Reff(ss),k_mealses);
            end
            os(jj)=cc; 
        end
    end
    Outbreak_County(ss,d_zero>0)=os;
end

end