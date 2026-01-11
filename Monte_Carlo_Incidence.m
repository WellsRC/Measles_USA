function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts,Imported_Case]=Monte_Carlo_Incidence(National_Annual_Reduction,NS,Scenario_Plot,Year_Reduced)
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat',"R_NHG","lambda_d",'k_mealses','lambda_RUCC');

[Imported_Case] = Case_Importation_Sample(Scenario_Plot,NS);

load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_i')


[~,RUCC_j] = Load_Transmission_Covariates(County_Data.GEOID);

clear County_Data


load(['National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated','Proportion_Age_Unvaccinated_Uninsured','Proportion_Age_Unvaccinated_Public','Proportion_Age_Unvaccinated_Private','Proportion_Age_Vaccinated_Uninsured','Proportion_Age_Vaccinated_Public','Proportion_Age_Vaccinated_Private');

Max_Outbreak=County_Data_Vaccine_Reduction.Total_Population(:).*(1-County_Data_Vaccine_Reduction.Total_Immunity(:));


Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

    
Case_Count(Reff<=1)=min(1./(1-Reff(Reff<1)),100);   
Case_Count(Reff>1 & Case_Count<=1+10^(-8))=1+10^(-8);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negatie hyper geomtric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_NHG=Max_Outbreak-1; % Substract one since we using the neg. hyprgeometric from o to maximal outbreak
N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcaultion of excess zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_0=zeros(size(Reff));
for cc=1:length(q_0)
    q_0(cc)=integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff(cc).*x)),0,1);
end

p_outbreak=(1-repmat(q_0,1,NS).^Imported_Case); % At least one of the imported cases triggers an utbreak
exp_case=repmat(Case_Count,1,NS).*p_outbreak; % Expected outbreak


p_RUCC=RUCC_j*lambda_RUCC;

z_ij=-lambda_d.*(Distance_Matrix_ij.^2);
w_ij=exp(z_ij); %1./(1+exp(-z_ij)); % Weight from population i (where the outbreak is) to population j
w_ij(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL

p_zero=zeros(size(Imported_Case));

parfor nn=1:size(p_zero,2)  
    
    Prev=repmat(exp_case(:,nn),1,length(exp_case(:,nn)))./Population_i;

    p_ij= exp(-repmat(p_RUCC',length(p_RUCC),1).*(1-repmat(q_0',size(w_ij,1),1)).*Prev.*w_ij); % Probability that county i does NOT trigger an outbeak in county j
    p_j = prod(p_ij,1)'; % Probability that an outbeak is NOT triggered in county j by domestic import
    p_zero(:,nn)=p_j.*(q_0.^Imported_Case(:,nn));
    
end

% Adjust the imported cases into Kansas to be discrete for the
% determinaition of vaccinated or unvaccianted cases. 


[Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,Proportion_Age_Unvaccinated_Uninsured,Proportion_Age_Unvaccinated_Public,Proportion_Age_Unvaccinated_Private,Proportion_Age_Vaccinated_Uninsured,Proportion_Age_Vaccinated_Public,Proportion_Age_Vaccinated_Private,NS,Imported_Case);

Contacts=repmat(County_Data_Vaccine_Reduction.All_Contacts,1,1,size(Vaccinated_Cases_County,3));
Total_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;

Contacts=repmat(County_Data_Vaccine_Reduction.Unvaccinated_Contacts,1,1,size(Vaccinated_Cases_County,3));
Unvaccinated_Contacts=Contacts.*(Vaccinated_Cases_County+Unvaccinated_Cases_County).*8;
end
