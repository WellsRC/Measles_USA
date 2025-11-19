function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts,Imported_Case]=Monte_Carlo_Incidence(National_Annual_Reduction,NS,Scenario_Plot,Year_Reduced)
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat','lambda_out',"R_NHG","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses');

[Imported_Case,Kansas_Discrete] = Case_Importation_Sample(Scenario_Plot,NS);

t_f= strcmp(County_Data.State,'Kansas');
Imported_Case(t_f,:)=Kansas_Discrete;
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i')

clear County_Data


load(['National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated');

Max_Outbreak=County_Data_Vaccine_Reduction.Total_Population(:).*(1-County_Data_Vaccine_Reduction.Total_Immunity(:));



[Case_Count,Reff]=Determine_Model_County_Case_Count(County_Data_Vaccine_Reduction,beta_j);

    
Case_Count(Reff<=1)=min(1./(1-Reff(Reff<1)),100);   
Case_Count(Reff>1 & Case_Count<=1+10^(-8))=1+10^(-8);

K_NHG=Max_Outbreak-1;
N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average

p_zero=zeros(size(Imported_Case));
for nn=1:size(p_zero,2)  
    [p_zero(:,nn),N_NHG,K_NHG]=Hurdle_Parameters(County_Data_Vaccine_Reduction,Case_Count,Max_Outbreak,R_NHG,Reff,k_mealses,Imported_Case(:,nn),lambda_0,lambda_i,lambda_j,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_out);
end
% Adjust the imported cases into Kansas to be discrete for the
% determinaition of vaccinated or unvaccianted cases. 


[Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,NS,Imported_Case);

Contacts=repmat(County_Data_Vaccine_Reduction.All_Contacts,1,1,size(Vaccinated_Cases_County,3));
Total_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;

Contacts=repmat(County_Data_Vaccine_Reduction.Unvaccinated_Contacts,1,1,size(Vaccinated_Cases_County,3));
Unvaccinated_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;
end
