function [p_zero,N_NHG,K_NHG]=Hurdle_Parameters(County_Data,Case_Count,Max_Outbreak,R_NHG,Reff,k_mealses,Imported_Case,lambda_0,lambda_i,lambda_j,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_out)
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

    p_outbreak=(1-q_0.^Imported_Case); % At least one of the imported cases triggers an utbreak
    exp_case=Case_Count(:).*p_outbreak(:); % Expected outbreak
    % Gravity model for flow from i to j
    z_ij=lambda_0+lambda_i.*log(Population_i)+lambda_j.*log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
    w_ij=1./(1+exp(-z_ij)); % Weight from population i (where the outbreak is) to population j
    w_ij(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL
    
    % https://pmc.ncbi.nlm.nih.gov/articles/PMC8521690/#sec21
    % We scale by the immunity level a a region with full immunity would
    % has greatest chance of n outbreak an with lowest immunuty the lowest
    % chance of no outbrak
    

    % Take the transpose of Total_Immunity a this is the desitantion where
    % the outbreak may possibely start as the outbreak is originating in i
    % and going to j
    p_ij= exp(-lambda_out.*(1-repmat(County_Data.Total_Immunity',size(w_ij,1),1)).*repmat(exp_case,1,size(w_ij,2)).*w_ij); % Probability that county i does NOT trigger an outbeak in county j

    p_j = prod(p_ij,1)'; % Probability that an outbeak is NOT triggered in county j by domestic import

    p_zero=p_j.*(q_0.^Imported_Case);

end