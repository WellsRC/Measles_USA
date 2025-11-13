function [Total_Cases_County]=Monte_Carlo_Incidence_REMOVE()
load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
NS=2500;
Scenario_Plot='Sample_2025';
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

load(['TEST_MODEL_IMM=60.mat'],'County_Data_Vaccine_Reduction');

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

[Total_Cases_County]=Monte_Carlo_Outbreak_County_REMOVE(F_NB,Max_Outbreak,p_c,Case_Count,k_nbin,Reff,k_mealses,NS,Imported_Case);

end
