function [p_zero,N_NHG,K_NHG]=Hurdle_Parameters(Case_Count,Max_Outbreak,R_NHG,Reff,Imported_Case,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_zero,k_mealses)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negatie hyper geomtric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K_NHG=Max_Outbreak-1; % Substract one since we using the neg. hyprgeometric from o to maximal outbreak
N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcaultion of excess zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NS=size(Imported_Case,2);

q_0=zeros(size(Reff));
for cc=1:length(q_0)
    q_0(cc)=integral(@(x)nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff(cc).*x)),0,1);
end

p_outbreak=(1-repmat(q_0,1,NS).^Imported_Case); % At least one of the imported cases triggers an utbreak
Prev=repmat(Case_Count./Population_i(:,1),1,NS).*p_outbreak; % Expected outbreak

z_ij=log(Population_i)+log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
w_ij=exp(z_ij); %1./(1+exp(-z_ij)); % Weight from population i (where the outbreak is) to population j
w_ij(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL
    

temp_zero=w_ij.*(1-repmat(q_0',size(w_ij,1),1));

p_j=exp(-lambda_zero.*(Prev')*temp_zero)';
p_zero=p_j.*(q_0.^Imported_Case);
end