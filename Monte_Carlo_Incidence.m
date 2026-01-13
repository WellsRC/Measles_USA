function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts,Imported_Case]=Monte_Carlo_Incidence(National_Annual_Reduction,NS,Scenario_Plot,Year_Reduced)
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

if(strcmp(Scenario_Plot,'Baseline'))
    load('Prior_log_Regression_Transmission.mat','prior_mean_transmission')
    [County_Transmission_X] = Load_Transmission_Covariates(County_Data.GEOID);
    
    [beta_j] = County_Transmission(prior_mean_transmission',County_Transmission_X.X);
    
    [~,Reff]=Determine_Model_County_Case_Count(County_Data,beta_j);
    
    t_f=strcmp(County_Data.County,'Gaines') & strcmp(County_Data.State,'Texas');
    
    q_0=integral(@(x)nbinpdf(0,0.23,0.23./(0.23+Reff(t_f).*x)),0,1);
    
    Import_Gaines=round(log(1-0.5)/log(q_0));    
    
else
    Import_Gaines=0;
end

[Imported_Case] = Case_Importation_Sample(Scenario_Plot,NS,Import_Gaines);

load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_i','Population_j')

load('Baseline_Estimate_Measles_Incidence.mat',"R_NHG","lambda_d",'k_mealses','lambda_zero');

clear County_Data


load(['National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated','Proportion_Age_Unvaccinated_Uninsured','Proportion_Age_Unvaccinated_Public','Proportion_Age_Unvaccinated_Private','Proportion_Age_Vaccinated_Uninsured','Proportion_Age_Vaccinated_Public','Proportion_Age_Vaccinated_Private');

Max_Outbreak=County_Data_Vaccine_Reduction.Total_Population(:).*(1-County_Data_Vaccine_Reduction.Total_Immunity(:));


Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

    
Case_Count(Reff<=1)=min(1./(1-Reff(Reff<1)),100);   
Case_Count(Reff>1 & Case_Count<=1+10^(-8))=1+10^(-8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcaultion of excess zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[p_zero,N_NHG,K_NHG]=Hurdle_Parameters(Case_Count,Max_Outbreak,R_NHG,Reff,Imported_Case,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_zero,k_mealses);

% Adjust the imported cases into Kansas to be discrete for the
% determinaition of vaccinated or unvaccianted cases. 



[Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,Proportion_Age_Unvaccinated_Uninsured,Proportion_Age_Unvaccinated_Public,Proportion_Age_Unvaccinated_Private,Proportion_Age_Vaccinated_Uninsured,Proportion_Age_Vaccinated_Public,Proportion_Age_Vaccinated_Private,NS,Imported_Case);

Contacts=repmat(County_Data_Vaccine_Reduction.All_Contacts,1,1,size(Vaccinated_Cases_County,3));
Total_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;

Contacts=repmat(County_Data_Vaccine_Reduction.Unvaccinated_Contacts,1,1,size(Vaccinated_Cases_County,3));
Unvaccinated_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;
end
