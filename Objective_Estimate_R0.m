function J = Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,X_Samp,L_Samp)

lambda_0=-10.^x(1);
lambda_i=10.^x(2);
lambda_j=10.^x(3);
lambda_d=10.^x(4);
k_mealses=10.^x(5); 
lambda_out=10.^x(6); 


indx_beta=round(x(7));
beta_j=Transmission_Relation(1-County_Data.Total_Immunity,X_Samp);

R_NHG=round(x(8));
Import_Gaines=round(x(9));
Import_Kansas=round(x(10));

L_Transmission=L_Samp;

K_NHG=Max_Outbreak;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source is unknown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume gains cunty as it was the epi center
t_f=strcmp(County_Data.County,'Gaines') & strcmp(County_Data.State,'Texas');
Imported_Case(t_f)=Imported_Case(t_f)+Import_Gaines; 

% Kanasa unrtain so distribute across state based on poulation
t_f= strcmp(County_Data.State,'Kansas');
Imported_Case(t_f)=Imported_Case(t_f)+Import_Kansas.*County_Data.Total_Population(t_f)./sum(County_Data.Total_Population(t_f)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Reproduction numbers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
Reff=zeros(length(beta_j),1);
for ss=1:length(beta_j)
    Reff(ss)=interp1(County_Data.beta_j',County_Data.R_eff(ss,:)',beta_j(ss))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Priors on Ro and Reeff and Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% https://static-content.springer.com/esm/art%3A10.1038%2Fnature04153/MediaObjects/41586_2005_BFnature04153_MOESM3_ESM.pdf
%  a=fmincon(@(x) norm(gaminv([0.05 0.95],x,0.23./(x-1))-[0.16 0.39]),11,[],[],[],[],1,10^3)
L_Measles=log(gampdf(k_mealses,11.5327,0.23/(11.5327-1)));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final sizes for the years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Case_Count=zeros(length(beta_j),1);

    for ss=1:length(beta_j)
        Case_Count(ss)=interp1(County_Data.beta_j',County_Data.Final_Size_Est(ss,:)',beta_j(ss))';
    end
    
    Case_Count(Reff<1)=1./(1-Reff(Reff<1));   
    
    Case_Count(Reff>=1 & Case_Count<=1+10^(-8))=1+10^(-8);


    N_NHG=(R_NHG.*K_NHG+(Case_Count-1).*K_NHG-(Case_Count-1))./(Case_Count-1); % Subtract one as a simplication to trunating to get the average
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    p_ij= exp(-lambda_out.*(1-repmat(County_Data.Total_Immunity,1,size(w_ij,2))).*repmat(exp_case,1,size(w_ij,2)).*w_ij); % Probability that county i does NOT trigger an outbeak in county j

    p_j = prod(p_ij,1)'; % Probability that an outbeak is NOT triggered in county j by domestic import
    
    
       
    p_zero=p_j.*(q_0.^Imported_Case);
    
    L_Known=zeros(size(p_zero));
    
    L_Unknown=zeros(size(p_zero));
    
    for cc=1:length(Known_Ind_Cases)
        if(isnan(Unknown_Ind_Cases_Weight(cc,1)))        
            if(Known_Ind_Cases(cc)==0)
                L_Known(cc)=log(p_zero(cc));
            elseif(Reff(cc)>1)
               L_Known(cc)=log((1-p_zero(cc)).*(1-neghyp_cdf(Known_Ind_Cases(cc)-1,N_NHG(cc),K_NHG(cc),R_NHG))) ; % as we decided to apprximate the trucnated distribution with a negative hypr geometric
            else
                temp_cdf=min(Chain_Size_Distribution_CDF(Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
                L_Known(cc)=log((1-p_zero(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen 
            end
        else
            if(isnan(Unknown_Ind_Cases(cc,2)))
                for uu=0:Unknown_Ind_Cases(cc,1)            
                    if(Known_Ind_Cases(cc)+uu==0)
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*p_zero(cc));
                    elseif(Reff(cc)>1)
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_zero(cc)).*(1-neghyp_cdf(Known_Ind_Cases(cc)+uu-1,N_NHG(cc),K_NHG(cc),R_NHG))); 
                    else
                        temp_cdf=min(Chain_Size_Distribution_CDF(uu+Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_zero(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
                    end
                end
            else
                for yy=0:Unknown_Ind_Cases(cc,2)            
                    for uu=0:Unknown_Ind_Cases(cc,1)            
                        if(Known_Ind_Cases(cc)+uu+yy==0)
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*p_zero(cc));
                        elseif(Reff(cc)>1)
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_zero(cc)).*(1-neghyp_cdf(Known_Ind_Cases(cc)+uu+yy-1,N_NHG(cc),K_NHG(cc),R_NHG))); 
                        else
                            temp_cdf=min(Chain_Size_Distribution_CDF(yy+uu+Known_Ind_Cases(cc)+1,Reff(cc),k_mealses),1);
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_zero(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
                        end
                    end
                end
            end
        end
    end
    
    if(~isinf(-mean(L_Unknown(:)+L_Known(:))-L_Measles -L_Transmission))
    
        [Outbreak_County_2025]=Monte_Carlo_Outbreak_County_Fitting(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,r_samp_pc_2025,r_samp_outbreak_2025);
        
        NOB_2025=sum(Outbreak_County_2025,1)'+sum(Imported_Case(:));    
        
        L_Weekly_2025=zeros(size(NOB_2025));
        x0_2025=log10([3.0129    6.9335    4.1784]);
        opts=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',5.*10^3,'FunctionTolerance',10^(-8),'StepTolerance',10^(-8),'MaxIterations',10^3);
        
        parfor jj=1:length(NOB_2025)
            [~,temp_L]=fmincon(@(g)-sum(log(nbinpdf(Week_Nat_Case_Count_2025(:)',10.^g(3),10.^g(3)./(10.^g(3)+(NOB_2025(jj)./gamcdf(52,10.^g(1),10.^g(2))).*(gamcdf(1:length(Week_Nat_Case_Count_2025),10.^g(1),10.^g(2))-gamcdf(0:(length(Week_Nat_Case_Count_2025)-1),10.^g(1),10.^g(2))))))),x0_2025,[],[],[],[],[-16 -16 -16],[3 3 3],[],opts);
            L_Weekly_2025(jj)=-temp_L;        
        end
    else
        L_Weekly_2025=0;
    end
    J=-mean(L_Unknown(:)+L_Known(:)) -mean(L_Weekly_2025(:)) -L_Measles -L_Transmission; 


end

