function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts,Imported_Case]=Monte_Carlo_Incidence(F_NB,National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6)
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat',"k_nbin","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses');


[Imported_Case] = Case_Importation_Sample(Scenario_Plot,NS);

load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i','County_GEOID')

indx_G=zeros(length(County_Data.County),1);
for cc=1:length(County_Data.County)
   indx_G(cc)=find(strcmp(County_GEOID,County_Data.GEOID{cc}));
end

Distance_Matrix_ij=Distance_Matrix_ij(indx_G,:);
Distance_Matrix_ij=Distance_Matrix_ij(:,indx_G);

Population_j=Population_j(indx_G,:);
Population_j=Population_j(:,indx_G);

Population_i=Population_i(indx_G,:);
Population_i=Population_i(:,indx_G);

FN_Age_Class={'Age_0_to_4','_Age_5_to_9','_Age_10_to_14','_Age_15_to_19','_Age_20_to_24'};
FN_Age_Class=FN_Age_Class(Age_Reduction);

if(~Age_0_to_6)
    load(['National_Reduction=' num2str(100*National_Reduction) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated');
else
    load(['National_Reduction=' num2str(100*National_Reduction) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated');
end
Reff_Seed=County_Data_Vaccine_Reduction.R_eff_Seed;
Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

Max_Outbreak=County_Data.Total_Population(:).*(1-County_Data.Total_Immunity(:));

Case_Count(Reff<1)=1./(1-Reff(Reff<1));

Reff_Seed=repmat(Reff_Seed,1,size(Imported_Case,2));


p_out=(1-nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^Imported_Case);
exp_case=repmat(Case_Count(:),1,NS).*p_out;
z_ij=lambda_0+lambda_i.*log(Population_i)+lambda_j.*log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
w_g=1./(1+exp(-z_ij')); % TOOK TRANSPOSE AS WE WANT THE FLOW WHERE THERE IS IMPORTED 
w_g(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL
Domestic_Import=w_g*exp_case;

p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^(Imported_Case+Domestic_Import);

[Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County(F_NB,Max_Outbreak,p_c,Case_Count,k_nbin,Reff,k_mealses,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,NS,Imported_Case);

Contacts=repmat(County_Data_Vaccine_Reduction.All_Contacts,1,1,size(Vaccinated_Cases_County,3));
Total_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;

Contacts=repmat(County_Data_Vaccine_Reduction.Unvaccinated_Contacts,1,1,size(Vaccinated_Cases_County,3));
Unvaccinated_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;
end
