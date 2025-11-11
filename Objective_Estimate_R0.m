function J = Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,F_NB,Max_Outbreak)
%,Importation_Cases_County_2023,Importation_Cases_County_2024,r_samp_pc_2023,r_samp_outbreak_2023,r_samp_pc_2024,r_samp_outbreak_2024,Week_Nat_Case_Count_2023,Week_Nat_Case_Count_2024)

beta_seed=10.^x(1);
lambda_0=x(2);
lambda_i=10.^x(3);
lambda_j=10.^x(4);
lambda_d=10.^x(5);
k_nbin=10.^x(6);
beta_j=10.^x(7);


% R0=interp1(County_Data.beta_j',County_Data.R_0',beta_seed)';
Reff_seed=interp1(County_Data.beta_j',County_Data.R_eff',beta_seed)';
Reff=interp1(County_Data.beta_j',County_Data.R_eff',beta_j)';


Case_Count=interp1(County_Data.beta_j',County_Data.Final_Size_Est',beta_seed)';

Case_Count(Reff<1)=1./(1-Reff(Reff<1));


k_mealses=0.23; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_out=(1-nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_seed)).^Imported_Case);
exp_case=Case_Count(:).*p_out(:);
z_ij=lambda_0+lambda_i.*log(Population_i)+lambda_j.*log(Population_j)-lambda_d.*log(Distance_Matrix_ij);
w_g=1./(1+exp(-z_ij')); % TOOK TRANSPOSE AS WE WANT THE FLOW WHERE THERE IS IMPORTED 
w_g(Distance_Matrix_ij==0)=0; % NO IMPACT ON DIAGONAL
Domestic_Import=w_g*exp_case;

p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_seed)).^(Imported_Case+Domestic_Import);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2023
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%     p_out=(1-nbinpdf(0,k_mealses,k_mealses./(k_mealses+repmat(Reff_Seed,1,size(Importation_Cases_County_2023,2)))).^Importation_Cases_County_2023);
%     exp_case=repmat(Case_Count(:),1,size(Importation_Cases_County_2023,2)).*p_out;
% 
%     Domestic_Import=w_g*exp_case;
%     p_c_2023=nbinpdf(0,k_mealses,k_mealses./(k_mealses+repmat(Reff_Seed,1,size(Importation_Cases_County_2023,2)))).^(Importation_Cases_County_2023+Domestic_Import);
% 
%     p_c_2023=mean(p_c_2023,2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 2024
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     p_out=(1-nbinpdf(0,k_mealses,k_mealses./(k_mealses+repmat(Reff_Seed,1,size(Importation_Cases_County_2023,2)))).^Importation_Cases_County_2024);
%     exp_case=repmat(Case_Count(:),1,size(Importation_Cases_County_2023,2)).*p_out;
% 
%     Domestic_Import=w_g*exp_case;
%     p_c_2024=nbinpdf(0,k_mealses,k_mealses./(k_mealses+repmat(Reff_Seed,1,size(Importation_Cases_County_2023,2)))).^(Importation_Cases_County_2024+Domestic_Import);
% 
%     p_c_2024=mean(p_c_2024,2);

L_Known=zeros(size(p_c));

L_Unknown=zeros(size(p_c));

for cc=1:length(Known_Ind_Cases)
    if(isnan(Unknown_Ind_Cases_Weight(cc,1)))        
        if(Known_Ind_Cases(cc)==0)
            L_Known(cc)=log(p_c(cc));
        elseif(Reff(cc)>1)
            z_nb=F_NB(log10(k_nbin),log10(max(Case_Count(cc),1.01)))';
            p_nb=1./(1+exp(-z_nb));
            pd = makedist('NegativeBinomial','R',k_nbin,'P',p_nb);
            pd = truncate(pd,1,inf);
            L_Known(cc)=log((1-p_c(cc)).*(1-cdf(pd,Known_Ind_Cases(cc)))); 
        else
            temp_cdf=min(Chain_Size_Distribution_CDF(Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
            L_Known(cc)=log((1-p_c(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen 
        end
    else
        if(isnan(Unknown_Ind_Cases(cc,2)))
            for uu=0:Unknown_Ind_Cases(cc,1)            
                if(Known_Ind_Cases(cc)+uu==0)
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*p_c(cc));
                elseif(Reff(cc)>1)
                    z_nb=F_NB(log10(k_nbin),log10(max(Case_Count(cc),1.01)))';
                    p_nb=1./(1+exp(-z_nb));
                    pd = makedist('NegativeBinomial','R',k_nbin,'P',p_nb);
                    pd = truncate(pd,1,inf);
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_c(cc)).*(1-cdf(pd,Known_Ind_Cases(cc)+uu))); 
                else
                    temp_cdf=min(Chain_Size_Distribution_CDF(uu+Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
                    L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_c(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
                end
            end
        else
            for yy=0:Unknown_Ind_Cases(cc,2)            
                for uu=0:Unknown_Ind_Cases(cc,1)            
                    if(Known_Ind_Cases(cc)+uu+yy==0)
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*p_c(cc));
                    elseif(Reff(cc)>1)
                        z_nb=F_NB(log10(k_nbin),log10(max(Case_Count(cc),1.01)))';
                        p_nb=1./(1+exp(-z_nb));
                        pd = makedist('NegativeBinomial','R',k_nbin,'P',p_nb);
                        pd = truncate(pd,1,inf);
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_c(cc)).*(1-cdf(pd,Known_Ind_Cases(cc)+yy+uu))); 
                    else
                        temp_cdf=min(Chain_Size_Distribution_CDF(yy+uu+Known_Ind_Cases(cc)+1,Reff(cc),k_mealses),1);
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_c(cc)).*(1-temp_cdf)); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
                    end
                end
            end
        end
    end
end
% R based on superpreading events https://jamanetwork.com/journals/jamapediatrics/fullarticle/2755836
% par=fmincon(@(x) norm(gaminv([0.025 0.5 0.975],x(1),x(2))-[5.0 6.1 18.1]),[5 2])
L_R0=log(gampdf(Reff(:),5.3110,1.6362)); 


% https://www.cdc.gov/mmwr/volumes/73/wr/mm7314a1.htm#:~:text=Among%20the%20338%20reported%20cases,rapidly%20investigating%20suspected%20measles%20cases.
% Table for transmission chains
% fmincon(@(x)-sum(log(betapdf([8 2 9 24 19]./[9 7 15 31 30],x(1),x(2)))),[3 1],[],[],[],[],[0 0],[5 5])
L_Control = log(betapdf(nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_seed)),3.2586,1.8919));
if(~isinf(-mean(L_Unknown(:)+L_Known(:))-mean(L_R0(:))))
    
    % [Outbreak_County_2023]=Monte_Carlo_Outbreak_County_Fitting(F_NB,Max_Outbreak,p_c_2023,Case_Count,k_nbin,Reff,k_mealses,r_samp_pc_2023,r_samp_outbreak_2023);
    % NOB_2023=sum(Outbreak_County_2023,1)'+mean(sum(Importation_Cases_County_2023,1));
    % 
    % [Outbreak_County_2024]=Monte_Carlo_Outbreak_County_Fitting(F_NB,Max_Outbreak,p_c_2024,Case_Count,k_nbin,Reff,k_mealses,r_samp_pc_2024,r_samp_outbreak_2024);
    % NOB_2024=sum(Outbreak_County_2024,1)'+mean(sum(Importation_Cases_County_2024,1));

    [Outbreak_County_2025]=Monte_Carlo_Outbreak_County_Fitting(F_NB,Max_Outbreak,p_c,Case_Count,k_nbin,Reff,k_mealses,r_samp_pc_2025,r_samp_outbreak_2025);
    NOB_2025=sum(Outbreak_County_2025,1)'+sum(Imported_Case(:));    
    
    L_Weekly_2025=zeros(size(NOB_2025));
    x0_2025=[3.0129    6.9335    4.1784];
    opts=optimoptions('fmincon','Display','none','MaxFunctionEvaluations',5.*10^3,'FunctionTolerance',10^(-8),'StepTolerance',10^(-8),'MaxIterations',10^3);

    parfor jj=1:length(NOB_2025)
        [~,temp_L]=fmincon(@(g)-sum(log(nbinpdf(Week_Nat_Case_Count_2025(:)',g(3),g(3)./(g(3)+(NOB_2025(jj)./gamcdf(52,g(1),g(2))).*(gamcdf(1:length(Week_Nat_Case_Count_2025),g(1),g(2))-gamcdf(0:(length(Week_Nat_Case_Count_2025)-1),g(1),g(2))))))),x0_2025,[],[],[],[],[0 0 0],[15 15 15],[],opts);
        L_Weekly_2025(jj)=-temp_L;        
    end
    % pd=fitdist(NOB_2024(:),'Kernel','Support','positive');
    % L_Weekly_2024=log(pdf(pd,sum(Week_Nat_Case_Count_2024)));
    % 
    % pd=fitdist(NOB_2023(:),'Kernel','Support','positive');
    % L_Weekly_2023=log(pdf(pd,sum(Week_Nat_Case_Count_2023)));
    % 
    % pd=fitdist(NOB_2025(:),'Kernel','Support','positive');
    % L_Weekly_2025=log(1-cdf(pd,sum(Week_Nat_Case_Count_2025)));
else
    L_Weekly_2025=0;
    % L_Weekly_2023=0;
    % L_Weekly_2024=0;
end
J=-mean(L_Unknown(:)+L_Known(:))-mean(L_R0(:))-mean(L_Weekly_2025(:)) -mean(L_Control(:)); % -L_Weekly_2023 -L_Weekly_2024;

end

