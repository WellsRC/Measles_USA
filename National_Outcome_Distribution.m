function [pd_cases,pd_hospital,pd_cost,pd_cost_per_case]=National_Outcome_Distribution(National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6)

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
[p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();
[Productivity_Cost_Under_15_Case,Productivity_Cost_15_plus_Case,Productivity_Cost_Under_15_Contact,Productivity_Cost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

Unvaccinated_Cases_County_Baseline=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
Vaccinated_Cases_County_Baseline=squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Productivity_loss_Cases=[Productivity_Cost_Under_15_Case.*ones(1,3) Productivity_Cost_15_plus_Case.*ones(1,15)];

Productivity_loss_Contacts=[Productivity_Cost_Under_15_Contact.*ones(1,3) Productivity_Cost_15_plus_Contact.*ones(1,15)];

Total_Productivity_loss_Cases=Productivity_loss_Cases*(Unvaccinated_Cases_County_Baseline+Vaccinated_Cases_County_Baseline);
Total_Productivity_loss_Contacts=Productivity_loss_Contacts*squeeze(sum(Unvaccinated_Contacts_Baseline,1));

Hospitalizations_Baseline=p_H_Unvaccinated*Unvaccinated_Cases_County_Baseline+p_H_Vaccinated*Vaccinated_Cases_County_Baseline;

Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,[1 2]));
% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Cost_Vaccination_Contacts=Cost_Vac_Age*squeeze(sum(Unvaccinated_Contacts_Baseline,1));
Cost_Case_Medical=Cost_per_Hospitalization.*Hospitalizations_Baseline+Cost_per_Non_Hospitalization.*(sum(Outbreak_Cases_County,1)-Hospitalizations_Baseline);
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;

Cost_Baseline=Contact_Tracing_Costs(:)+Cost_Case_Medical(:)+Testing_Cost(:)+Cost_Vaccination_Contacts(:)+Total_Productivity_loss_Cases(:)+Total_Productivity_loss_Contacts(:);

% Cases    
temp_c=sum(Outbreak_Cases_County,1);
pd_cases=fitdist(temp_c(:),'Kernel','Support','positive');

% Hospital
temp_c=Hospitalizations_Baseline;
pd_hospital=fitdist(temp_c(:),'Kernel','Support','positive');

% Cost
temp_c=Cost_Baseline./10^6;
pd_cost=fitdist(temp_c(:),'Kernel','Support','positive');

% Cost per case
Cost_per_Case=Cost_Baseline./(sum(Outbreak_Cases_County,1)');
temp_c=Cost_per_Case./10^3;
pd_cost_per_case=fitdist(temp_c(:),'Kernel','Support','positive');

end
