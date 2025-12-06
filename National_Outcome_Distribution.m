function [pd_cases,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_med_cost_uninsured,pd_med_cost_public,pd_med_cost_private,pd_test_vac_cost,pd_ct_cost,pd_outbreak_response_cost]=National_Outcome_Distribution(National_Annual_Reduction,Scenario_Plot,Year_Reduced)

[p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization]=Hospitalization_Probability();
[Productivity_Days_Lost_Under_15_Case,Productivity_Days_Lost_15_plus_Case,Productivity_Days_Lost_Under_15_Contact,Productivity_Days_Lost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

load(['National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction')
    

[Cost_Case_Medical_Uninsured]= Compute_Direct_Medical_Costs('Uninsured',County_Data_Vaccine_Reduction,Scenario_Plot,National_Annual_Reduction,Year_Reduced,p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization,Cost_per_Non_Hospitalization);
[Cost_Case_Medical_Public]= Compute_Direct_Medical_Costs('Public',County_Data_Vaccine_Reduction,Scenario_Plot,National_Annual_Reduction,Year_Reduced,p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization,Cost_per_Non_Hospitalization);
[Cost_Case_Medical_Private]= Compute_Direct_Medical_Costs('Private',County_Data_Vaccine_Reduction,Scenario_Plot,National_Annual_Reduction,Year_Reduced,p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization,Cost_per_Non_Hospitalization);

[Total_Productivity_loss_Cases,Total_Productivity_loss_Contacts]=Compute_Productivity_Losses(County_Data_Vaccine_Reduction,Scenario_Plot,National_Annual_Reduction,Year_Reduced,Productivity_Days_Lost_Under_15_Case,Productivity_Days_Lost_15_plus_Case,Productivity_Days_Lost_Under_15_Contact,Productivity_Days_Lost_15_plus_Contact);
load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline');

Hospitalizations_Baseline=p_H_Unvaccinated*squeeze(sum(Unvaccinated_Cases_County_Baseline,1))+p_H_Vaccinated*squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,[1 2]));
% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Cost_Vaccination_Contacts=Cost_Vac_Age*squeeze(sum(Unvaccinated_Contacts_Baseline,1));

Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;


Total_Productivity_loss=Total_Productivity_loss_Cases(:)+Total_Productivity_loss_Contacts(:);
Direct_Medical_Costs=Cost_Case_Medical_Uninsured(:)+Cost_Case_Medical_Public(:)+Cost_Case_Medical_Private(:);
Testing_Vaccination_Contacts_Cost=Testing_Cost(:)+Cost_Vaccination_Contacts(:);
Contact_Tracing_Costs=Contact_Tracing_Costs(:);



Cost_Baseline=Contact_Tracing_Costs(:)+Direct_Medical_Costs(:)+Testing_Cost(:)+Cost_Vaccination_Contacts(:)+Total_Productivity_loss_Cases(:)+Total_Productivity_loss_Contacts(:);

temp_c=Contact_Tracing_Costs(:)+Testing_Cost(:)+Cost_Vaccination_Contacts(:);
pd_outbreak_response_cost=fitdist(temp_c(:),'Kernel','Support','positive');


% Cases    
temp_c=sum(Total_Cases_County,1);
pd_cases=fitdist(temp_c(:),'Kernel','Support','positive');

% Hospital
temp_c=Hospitalizations_Baseline;
pd_hospital=fitdist(temp_c(:),'Kernel','Support','positive');

% Cost
temp_c=Cost_Baseline./10^6;
pd_cost=fitdist(temp_c(:),'Kernel','Support','positive');

% Cost per case
Cost_per_Case=Cost_Baseline./(sum(Total_Cases_County,1)');
temp_c=Cost_per_Case./10^3;
pd_cost_per_case=fitdist(temp_c(:),'Kernel','Support','positive');

% Productivity_loss    
temp_c=Total_Productivity_loss;
pd_pro_loss=fitdist(temp_c(:),'Kernel','Support','positive');

% Medical_Costs
temp_c=Direct_Medical_Costs;
pd_med_cost=fitdist(temp_c(:),'Kernel','Support','positive');

% Medical_Costs
temp_c=Cost_Case_Medical_Uninsured;
temp_c(temp_c==0)=10^(-16);
pd_med_cost_uninsured=fitdist(temp_c(:),'Kernel','Support','positive');

% Medical_Costs
temp_c=Cost_Case_Medical_Public;
temp_c(temp_c==0)=10^(-16);
pd_med_cost_public=fitdist(temp_c(:),'Kernel','Support','positive');

% Medical_Costs
temp_c=Cost_Case_Medical_Private;
temp_c(temp_c==0)=10^(-16);
pd_med_cost_private=fitdist(temp_c(:),'Kernel','Support','positive');

% Testing_Vaccination
temp_c=Testing_Vaccination_Contacts_Cost;
pd_test_vac_cost=fitdist(temp_c(:),'Kernel','Support','positive');

% Contact_Tracing_Costs
temp_c=Contact_Tracing_Costs;
pd_ct_cost=fitdist(temp_c(:),'Kernel','Support','positive');


end
