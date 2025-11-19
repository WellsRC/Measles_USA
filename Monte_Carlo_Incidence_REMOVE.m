function [Total_Cases_County]=Monte_Carlo_Incidence_REMOVE()

NS=500;
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat','lambda_out',"R_NHG","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses');

[Imported_Case,Kansas_Discrete] = Case_Importation_Sample('Baseline',NS);

t_f= strcmp(County_Data.State,'Kansas');
Imported_Case(t_f,:)=Kansas_Discrete;
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i')

clear County_Data

load('TEST_MODEL_IMM=60.mat');

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


[Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,NS,Imported_Case);

end
