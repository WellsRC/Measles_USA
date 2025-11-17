function [Total_Cases_County]=Monte_Carlo_Incidence_REMOVE()

NS=500;
Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat','lambda_out',"R_NHG","lambda_0","lambda_i","lambda_j","lambda_d",'k_mealses');

[Imported_Case,Kansas_Discrete] = Case_Importation_Sample('Baseline',NS);

t_f= strcmp(County_Data.State,'Kansas');
Imported_Case(t_f,:)=Kansas_Discrete;
load('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i','County_GEOID')

indx_G=zeros(length(County_Data.County),1);
for cc=1:length(County_Data.County)
   indx_G(cc)=find(strcmp(County_GEOID,County_Data.GEOID{cc}));
end

clear County_Data
Distance_Matrix_ij=Distance_Matrix_ij(indx_G,:);
Distance_Matrix_ij=Distance_Matrix_ij(:,indx_G);

Population_j=Population_j(indx_G,:);
Population_j=Population_j(:,indx_G);

Population_i=Population_i(indx_G,:);
Population_i=Population_i(:,indx_G);

load('TEST_MODEL_IMM=60.mat');

Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;
Max_Outbreak=County_Data_Vaccine_Reduction.Total_Population(:).*(1-County_Data_Vaccine_Reduction.Total_Immunity(:));
Case_Count(Reff<1)=1./(1-Reff(Reff<1));
Case_Count(Reff>=1 & Case_Count<=1+10^(-8))=1+10^(-8);

K_NHG=Max_Outbreak;
N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average

f_correct=N_NHG-R_NHG-K_NHG<0;
N_NHG(f_correct)=N_NHG(f_correct)+abs(min(N_NHG-R_NHG-K_NHG));
p_zero=zeros(size(Imported_Case));

q_0=zeros(size(Reff));
for cc=1:length(q_0)
    q_0(cc)=integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff(cc).*x)),0,1);
end

% Gravity model for flow from i to j
z_ij=lambda_0+lambda_i.*log(Population_i)+lambda_j.*log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
w_ij=1./(1+exp(-z_ij)); % Weight from population i (where the outbreak is) to population j
w_ij(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL
for nn=1:size(p_zero,2)  
    p_outbreak=(1-q_0.^Imported_Case(:,nn)); % At least one of the imported cases triggers an utbreak
    exp_case=Case_Count(:).*p_outbreak(:); % Expected outbreak
    
    
    % https://pmc.ncbi.nlm.nih.gov/articles/PMC8521690/#sec21
    % We scale by the immunity level a a region with full immunity would
    % has greatest chance of n outbreak an with lowest immunuty the lowest
    % chance of no outbrak
    
    
    % Take the transpose of Total_Immunity a this is the desitantion where
    % the outbreak may possibely start as the outbreak is originating in i
    % and going to j
    p_ij= exp(-lambda_out.*(1-repmat(County_Data_Vaccine_Reduction.Total_Immunity',size(w_ij,1),1)).*repmat(exp_case,1,size(w_ij,2)).*w_ij); % Probability that county i does NOT trigger an outbeak in county j
    
    p_j = prod(p_ij,1)'; % Probability that an outbeak is NOT triggered in county j by domestic import
    
    p_zero(:,nn)=p_j.*(q_0.^Imported_Case(:,nn));
end
% Adjust the imported cases into Kansas to be discrete for the
% determinaition of vaccinated or unvaccianted cases. 


[Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,NS,Imported_Case);

end
