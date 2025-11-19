function J = Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Week_Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,X_Samp,L_Samp)

lambda_0=-10.^x(1);
lambda_i=10.^x(2);
lambda_j=10.^x(3);
lambda_d=10.^x(4);
k_mealses=10.^x(5); 
lambda_out=10.^x(6); 


indx_beta=round(x(7));
beta_j=Transmission_Relation(1-County_Data.Total_Immunity,X_Samp(indx_beta,:));

R_NHG=round(x(8));
Import_Gaines=round(x(9));
Import_Kansas=round(x(10));

L_Transmission=L_Samp(indx_beta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source is unknown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume gains cunty as it was the epi center
t_f=strcmp(County_Data.County,'Gaines') & strcmp(County_Data.State,'Texas');
Imported_Case(t_f)=Imported_Case(t_f)+Import_Gaines; 

% Kanasa unrtain so distribute across state based on poulation
t_f= strcmp(County_Data.State,'Kansas');
Imported_Case(t_f)=Imported_Case(t_f)+Import_Kansas.*County_Data.Total_Population(t_f)./sum(County_Data.Total_Population(t_f)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Priors on Ro and Reeff and Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% https://static-content.springer.com/esm/art%3A10.1038%2Fnature04153/MediaObjects/41586_2005_BFnature04153_MOESM3_ESM.pdf
%  a=fmincon(@(x) norm(gaminv([0.05 0.95],x,0.23./(x-1))-[0.16 0.39]),11,[],[],[],[],1,10^3)
L_Measles=log(gampdf(k_mealses,11.5327,0.23/(11.5327-1)));

[Case_Count,Reff]=Determine_Model_County_Case_Count(County_Data,beta_j);
[p_zero,N_NHG,K_NHG]=Hurdle_Parameters(County_Data,Case_Count,Max_Outbreak,R_NHG,Reff,k_mealses,Imported_Case,lambda_0,lambda_i,lambda_j,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_out);
    
    L_Known=zeros(size(p_zero));
    
    L_Unknown=zeros(size(p_zero));
    
    for cc=1:length(Known_Ind_Cases)
        if(isnan(Unknown_Ind_Cases_Weight(cc,1)))        
            if(Known_Ind_Cases(cc)==0)
                L_Known(cc)=log(p_zero(cc));
            elseif(Reff(cc)>1)
               L_Known(cc)=log((1-p_zero(cc)).*(neghyp_pdf(Known_Ind_Cases(cc)-1,N_NHG(cc),K_NHG(cc),R_NHG))) ; % as we decided to apprximate the trucnated distribution with a negative hypr geometric
            else
                % temp_cdf=min(Chain_Size_Distribution_CDF(Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
                L_Known(cc)=log((1-p_zero(cc)).*(Chain_Size_Distribution(Known_Ind_Cases(cc),Reff(cc),k_mealses))); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen 
            end
        else
            if(isnan(Unknown_Ind_Cases(cc,2)))
                for uu=0:Unknown_Ind_Cases(cc,1)            
                    if(Known_Ind_Cases(cc)+uu==0)
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*p_zero(cc));
                    elseif(Reff(cc)>1)
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_zero(cc)).*(neghyp_pdf(Known_Ind_Cases(cc)+uu-1,N_NHG(cc),K_NHG(cc),R_NHG))); 
                    else
                        % temp_cdf=min(Chain_Size_Distribution_CDF(uu+Known_Ind_Cases(cc),Reff(cc),k_mealses),1);
                        L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-p_zero(cc)).*(Chain_Size_Distribution(uu+Known_Ind_Cases(cc),Reff(cc),k_mealses))); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
                    end
                end
            else
                for yy=0:Unknown_Ind_Cases(cc,2)            
                    for uu=0:Unknown_Ind_Cases(cc,1)            
                        if(Known_Ind_Cases(cc)+uu+yy==0)
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*p_zero(cc));
                        elseif(Reff(cc)>1)
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_zero(cc)).*(neghyp_pdf(Known_Ind_Cases(cc)+uu+yy-1,N_NHG(cc),K_NHG(cc),R_NHG))); 
                        else
                            % temp_cdf=min(Chain_Size_Distribution_CDF(yy+uu+Known_Ind_Cases(cc)+1,Reff(cc),k_mealses),1);
                            L_Unknown(cc)=L_Unknown(cc)+log(binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-p_zero(cc)).*(Chain_Size_Distribution(yy+uu+Known_Ind_Cases(cc)+1,Reff(cc),k_mealses))); %add one to known cases as assuming there was an introduction from some place for this local transmission to happen
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
            [~,temp_L]=fmincon(@(g)-sum(log(nbinpdf(Week_Nat_Case_Count_2025(:)',10.^g(3),10.^g(3)./(10.^g(3)+(NOB_2025(jj)./gamcdf(52,10.^g(1),10.^g(2))).*(gamcdf(1:length(Week_Nat_Case_Count_2025),10.^g(1),10.^g(2))-gamcdf(0:(length(Week_Nat_Case_Count_2025)-1),10.^g(1),10.^g(2))))))),x0_2025,[],[],[],[],[-16 -16 log10(5)],[3 3 3],[],opts);
            L_Weekly_2025(jj)=-temp_L;        
        end
    else
        L_Weekly_2025=0;
    end
    J=-mean(L_Unknown(:)+L_Known(:)) -mean(L_Weekly_2025(:)) -L_Measles -L_Transmission;


end

