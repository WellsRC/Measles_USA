function J = Objective_Estimate_R0(x,County_Data,Imported_Case,Known_Ind_Cases,Unknown_Ind_Cases,Unknown_Ind_Cases_Weight,Population_i,Population_j,Distance_Matrix_ij,Nat_Case_Count_2025,r_samp_pc_2025,r_samp_outbreak_2025,Max_Outbreak,County_Transmission_X)

lambda_d=10.^x(1);
k_mealses=10.^x(2); 

R_NHG=round(x(3));

beta_transmission=x(4:32);

lambda_zero=10.^x(33);

[beta_j] = County_Transmission(beta_transmission,County_Transmission_X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Priors on Ro and Reeff and Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% https://static-content.springer.com/esm/art%3A10.1038%2Fnature04153/MediaObjects/41586_2005_BFnature04153_MOESM3_ESM.pdf
%  a=fmincon(@(x) norm(gaminv([0.05 0.95],x,0.23./(x-1))-[0.16 0.39]),11,[],[],[],[],1,10^3)
L_Measles=log(gampdf(k_mealses,11.5327,0.23/(11.5327-1)));

[Case_Count,Reff]=Determine_Model_County_Case_Count(County_Data,beta_j);
[p_zero,N_NHG,K_NHG]=Hurdle_Parameters(Case_Count,Max_Outbreak,R_NHG,Reff,Imported_Case,Population_i,Population_j,lambda_d,Distance_Matrix_ij,lambda_zero,k_mealses);
avg_p_zero=mean(p_zero,2);
L_Known=zeros(size(avg_p_zero));
    
    L_Unknown=zeros(size(avg_p_zero));
    
    for cc=1:length(Known_Ind_Cases)
        if(isnan(Unknown_Ind_Cases_Weight(cc,1)))        
            if(Known_Ind_Cases(cc)==0)
                L_Known(cc)=log(avg_p_zero(cc));
            elseif(Reff(cc)>1)
               L_Known(cc)=log((1-avg_p_zero(cc)))+log((neghyp_pdf(Known_Ind_Cases(cc)-1,N_NHG(cc),K_NHG(cc),R_NHG))) ; % as we decided to apprximate the trucnated distribution with a negative hypr geometric
            else
                L_Known(cc)=log((1-avg_p_zero(cc)))+log((Chain_Size_Distribution(Known_Ind_Cases(cc),Reff(cc),k_mealses))); 
            end
        else
            if(isnan(Unknown_Ind_Cases(cc,2)))
                for uu=0:Unknown_Ind_Cases(cc,1)            
                    if(Known_Ind_Cases(cc)+uu==0)
                        L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*avg_p_zero(cc);
                    elseif(Reff(cc)>1)
                        L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-avg_p_zero(cc)).*neghyp_pdf(Known_Ind_Cases(cc)+uu-1,N_NHG(cc),K_NHG(cc),R_NHG); 
                    else
                        L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*(1-avg_p_zero(cc)).*(Chain_Size_Distribution(uu+Known_Ind_Cases(cc),Reff(cc),k_mealses));
                    end
                end
            else
                for yy=0:Unknown_Ind_Cases(cc,2)            
                    for uu=0:Unknown_Ind_Cases(cc,1)            
                        if(Known_Ind_Cases(cc)+uu+yy==0)
                            L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*avg_p_zero(cc);
                        elseif(Reff(cc)>1)
                            L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-avg_p_zero(cc)).*neghyp_pdf(Known_Ind_Cases(cc)+uu+yy-1,N_NHG(cc),K_NHG(cc),R_NHG); 
                        else
                            L_Unknown(cc)=L_Unknown(cc)+binopdf(uu,Unknown_Ind_Cases(cc,1),Unknown_Ind_Cases_Weight(cc,1)).*binopdf(yy,Unknown_Ind_Cases(cc,2),Unknown_Ind_Cases_Weight(cc,2)).*(1-avg_p_zero(cc)).*Chain_Size_Distribution(yy+uu+Known_Ind_Cases(cc),Reff(cc),k_mealses); 
                        end
                    end
                end
            end
        end
    end
    L_Unknown=log(L_Unknown);
    L_Unknown(isnan(Unknown_Ind_Cases_Weight(:,1)))=0;

    if(~isinf(-mean(L_Unknown(:)+L_Known(:))-L_Measles) && ~isnan(-mean(L_Unknown(:)+L_Known(:))-L_Measles))
        try 
            [Outbreak_County_2025]=Monte_Carlo_Outbreak_County_Fitting(Max_Outbreak,p_zero,N_NHG,K_NHG,R_NHG,Reff,k_mealses,r_samp_pc_2025,r_samp_outbreak_2025);
    
            NOB_2025=sum(Outbreak_County_2025+Imported_Case,1)';    
    

            min_Case=unique(sum(Imported_Case,1));
            NOB_2025(NOB_2025==min_Case)=min_Case+10^(-8);
            pdf_nat=fitdist(NOB_2025(:),'Kernel','Support',[min_Case 10^9]);

            L_Nat_2025=log(pdf(pdf_nat,Nat_Case_Count_2025)); % Want the mode of the dist to be scaled equal for each distribution as we will use that for the central measure
        catch
            L_Nat_2025=NaN;
        end
    else
        L_Nat_2025=0;
    end


    % temp_Case=Known_Ind_Cases;
    % temp_Case(isnan(temp_Case))=0;
    % 
    % Affected_Counties=temp_Case;
    % Affected_Counties(Affected_Counties>1)=1;
    % 
    % Affected_Counties=sum(Affected_Counties);

    % Number_Outbreaks=49;
    J=-sum(L_Unknown(:)+L_Known(:))./sum(Known_Ind_Cases>=3) -L_Measles -L_Nat_2025; % Scaled county likelihood such that it did not outweight national incidence

    if(~isreal(J))
        J=NaN;
    end

end

