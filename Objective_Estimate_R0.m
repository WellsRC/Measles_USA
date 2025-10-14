function J = Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,G,Nat_Case_Count,epi_weight)

beta_seed=10.^x(1);
lambda_g=-10.^x(2);
k_nbin=10.^x(3);
beta_j=10.^x(4);

R0=interp1(County_Data.beta_j',County_Data.R_0',beta_seed)';
Reff_Seed=interp1(County_Data.beta_j',County_Data.R_eff',beta_seed)';
Reff=interp1(County_Data.beta_j',County_Data.R_eff',beta_j)';

Case_Count=interp1(County_Data.beta_j',County_Data.Final_Size_Est',beta_j)';

Case_Count(Reff<1)=1./(1-Reff(Reff<1));

GM=G*Case_Count(:)./(length(Reff(:))-1);

k_mealses=0.23; 
p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^(Imported_Case+exp(lambda_g.*GM));

L_Known=zeros(size(p_c));

L_Unknown=zeros(size(p_c));

for cc=1:length(Known_Ind_Cases)
    if(isnan(Unknown_Ind_Cases_Weight(cc)))
        if(Reff(cc)>1)
            if(Known_Ind_Cases(cc)==0)
                L_Known(cc)=log(p_c(cc)+(1-p_c(cc)).*nbinpdf(Known_Ind_Cases(cc),k_nbin,k_nbin./(k_nbin+Case_Count(cc).*epi_weight(cc)))); %log(p_c(cc));%
            else
                L_Known(cc)=log((1-p_c(cc)).*nbinpdf(Known_Ind_Cases(cc),k_nbin,k_nbin./(k_nbin+Case_Count(cc).*epi_weight(cc)))); %log(1-p_c(cc));%
            end
        else
            if(Known_Ind_Cases(cc)==0)
                L_Known(cc)=log(p_c(cc)+(1-p_c(cc)).*Chain_Size_Distribution(Known_Ind_Cases(cc)+1,Reff(cc),k_mealses)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen (i.e. for zero cases there was only the introution and nothing else) %log(p_c(cc));%
            else
                L_Known(cc)=log((1-p_c(cc)).*Chain_Size_Distribution(Known_Ind_Cases(cc),Reff(cc),k_mealses)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen %log(1-p_c(cc));%
            end
        end
    else
        for uu=0:Unknown_Ind_Cases(cc)
            if(Reff(cc)>1)
                if(Known_Ind_Cases(cc)+uu==0)
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(p_c(cc)+(1-p_c(cc)).*nbinpdf(Known_Ind_Cases(cc)+uu,k_nbin,k_nbin./(k_nbin+Case_Count(cc).*epi_weight(cc))))); %L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*p_c(cc)); %
                else
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(1-p_c(cc)).*nbinpdf(Known_Ind_Cases(cc)+uu,k_nbin,k_nbin./(k_nbin+Case_Count(cc).*epi_weight(cc)))); %L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(1-p_c(cc))); %
                end
            else
                if(Known_Ind_Cases(cc)+uu==0)
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(p_c(cc)+(1-p_c(cc)).*Chain_Size_Distribution(1+Known_Ind_Cases(cc)+uu,Reff(cc),k_mealses))); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen (i.e. for zero cases there was only the introution and nothing else) %L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*p_c(cc)); %
                else
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(1-p_c(cc)).*Chain_Size_Distribution(1+Known_Ind_Cases(cc)+uu,Reff(cc),k_mealses)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen  %L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc),Unknown_Ind_Cases_Weight(cc)).*(1-p_c(cc))); %
                end
            end
        end
    end
end

[Outbreak_County]=Monte_Carlo_Outbreak_County(p_c,Case_Count,k_nbin,Reff,k_mealses,250);
NOB=sum(Outbreak_County,1)';
NOB(NOB==0)=10^(-16);

logn_p=lognfit(NOB);
L_Nat=log(lognpdf(Nat_Case_Count,logn_p(1),logn_p(2)));


% L_Nat(end)=log(1-logncdf(Nat_Case_Count(end),logn_p(1),logn_p(2))); % the last point is censored as the year is not over
% R based on superpreading events https://jamanetwork.com/journals/jamapediatrics/fullarticle/2755836
% par=fmincon(@(x) norm(gaminv([0.025 0.5 0.975],x(1),x(2))-[5.0 6.1 18.1]),[5 2])
L_R0=gampdf(R0(:),5.3110,1.6362); 

J=-mean(L_Unknown(:)+L_Known(:))-sum(L_Nat(:))-mean(L_R0(:));

end

