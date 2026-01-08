function [Case_Count,Reff,R0]=Determine_Model_County_Case_Count(County_Data,beta_j)

Reff=zeros(size(County_Data.R_eff,1),1);
R0=zeros(size(County_Data.R_eff,1),1);
Case_Count=zeros(size(County_Data.R_eff,1),1);

for ss=1:size(County_Data.R_eff,1)
    Reff(ss)=interp1(County_Data.beta_j',County_Data.R_eff(ss,:)',beta_j(ss))';
    R0(ss)=interp1(County_Data.beta_j',County_Data.R_0(ss,:)',beta_j(ss))';
    Case_Count(ss)=interp1(County_Data.beta_j',County_Data.Final_Size_Est(ss,:)',beta_j(ss))';
    if(isnan(Reff(ss)))

        Reff(ss)=pchip(County_Data.beta_j',County_Data.R_eff(ss,:)',beta_j(ss))';
        R0(ss)=pchip(County_Data.beta_j',County_Data.R_0(ss,:)',beta_j(ss))';
        t_f=find(County_Data.Final_Size_Est(ss,:)'>0);
        if(~isempty(t_f))
            Case_Count(ss)=pchip(County_Data.beta_j(t_f)',County_Data.Final_Size_Est(ss,t_f)',beta_j(ss))';
        else
            Case_Count(ss)=0;
        end
    end
end

Case_Count(Reff<=1)=min(1./(1-Reff(Reff<1)),100);   

Case_Count(Reff>1 & Case_Count<=1+10^(-8))=1+10^(-8);
end