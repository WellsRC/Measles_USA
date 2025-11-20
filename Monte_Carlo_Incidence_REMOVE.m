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


Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

    
Case_Count(Reff<=1)=min(1./(1-Reff(Reff<1)),100);   
Case_Count(Reff>1 & Case_Count<=1+10^(-8))=1+10^(-8);


K_NHG=Max_Outbreak-1;
N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average


q_0=zeros(size(Reff));
for cc=1:length(q_0)
    q_0(cc)=integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff(cc).*x)),0,1);
end

p_outbreak=(1-repmat(q_0,1,NS).^Imported_Case); % At least one of the imported cases triggers an utbreak
exp_case=repmat(Case_Count,1,NS).*p_outbreak; % Expected outbreak
% Gravity model for flow from i to j
z_ij=lambda_0+lambda_i.*log(Population_i)+lambda_j.*log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
w_ij=1./(1+exp(-z_ij)); % Weight from population i (where the outbreak is) to population j
w_ij(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL

p_zero=zeros(size(Imported_Case));

Immunity=County_Data_Vaccine_Reduction.Total_Immunity;
parfor nn=1:size(p_zero,2)  

    % https://pmc.ncbi.nlm.nih.gov/articles/PMC8521690/#sec21
% We scale by the immunity level a a region with full immunity would
% has greatest chance of n outbreak an with lowest immunuty the lowest
% chance of no outbrak

    % Take the transpose of Total_Immunity a this is the desitantion where
    % the outbreak may possibely start as the outbreak is originating in i
    % and going to j
    p_ij= exp(-lambda_out.*(1-repmat(Immunity',size(w_ij,1),1)).*repmat(exp_case(:,nn),1,size(w_ij,2)).*w_ij); % Probability that county i does NOT trigger an outbeak in county j
    p_j = prod(p_ij,1)'; % Probability that an outbeak is NOT triggered in county j by domestic import
    p_zero(:,nn)=p_j.*(q_0.^Imported_Case(:,nn));
    
end





[Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,NS,Imported_Case);

end
