function [Cost_per_Case_Difference]=Compute_Change_Cost_per_Case(Age_Reduction,NS,Scenario_Plot,Age_0_to_6)

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
[p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();
[Productivity_Cost_Under_15_Case,Productivity_Cost_15_plus_Case,Productivity_Cost_Under_15_Contact,Productivity_Cost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,0,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

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

Cost_per_Case_Baseline=Cost_Baseline./sum(Outbreak_Cases_County',2);


[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,0.01,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

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

Cost_per_Case_Reduction_1=Cost_Baseline./sum(Outbreak_Cases_County',2);


[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,0.025,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

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

Cost_per_Case_Reduction_2_5=Cost_Baseline./sum(Outbreak_Cases_County',2);

[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,0.05,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

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

Cost_per_Case_Reduction_5=Cost_Baseline./sum(Outbreak_Cases_County',2);

DX_1=Cost_per_Case_Reduction_1-Cost_per_Case_Baseline;
DX_2_5=Cost_per_Case_Reduction_2_5-Cost_per_Case_Baseline;
DX_5=Cost_per_Case_Reduction_5-Cost_per_Case_Baseline;

pd_cost_per_case_reduction_1=fitdist(DX_1(:),'Kernel');
pd_cost_per_case_reduction_2_5=fitdist(DX_2_5(:),'Kernel');
pd_cost_per_case_reduction_5=fitdist(DX_5(:),'Kernel');


Scenario={'1% Reduction'; '2.5% Reduction'; '5% Reduction'};
Cost_Difference=cell(3,1);
Cost_Difference_RAW=cell(3,1);
% Costs
    x0=icdf(pd_cost_per_case_reduction_1,0.5);
    lb=icdf(pd_cost_per_case_reduction_1,0.05);
    ub=icdf(pd_cost_per_case_reduction_1,0.95);
    temp_X=fmincon(@(z)-log(pdf(pd_cost_per_case_reduction_1,z)),(x0),[],[],[],[],(lb),(ub));
    Cost_Difference{1} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_cost_per_case_reduction_1,0.025),'%4.2f') char(8211) num2str(icdf(pd_cost_per_case_reduction_1,0.975),'%4.2f') ')'];
    Cost_Difference_RAW{1} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(prctile(DX_1(:),2.5),'%4.2f') char(8211) num2str(prctile(DX_1(:),97.5),'%4.2f') ')'];

    x0=icdf(pd_cost_per_case_reduction_2_5,0.5);
    lb=icdf(pd_cost_per_case_reduction_2_5,0.05);
    ub=icdf(pd_cost_per_case_reduction_2_5,0.95);
    temp_X=fmincon(@(z)-log(pdf(pd_cost_per_case_reduction_2_5,z)),(x0),[],[],[],[],(lb),(ub));
    Cost_Difference{2} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_cost_per_case_reduction_2_5,0.025),'%4.2f') char(8211) num2str(icdf(pd_cost_per_case_reduction_2_5,0.975),'%4.2f') ')'];
    Cost_Difference_RAW{2} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(prctile(DX_2_5(:),2.5),'%4.2f') char(8211) num2str(prctile(DX_2_5(:),97.5),'%4.2f') ')'];

    x0=icdf(pd_cost_per_case_reduction_5,0.5);
    lb=icdf(pd_cost_per_case_reduction_5,0.05);
    ub=icdf(pd_cost_per_case_reduction_5,0.95);
    temp_X=fmincon(@(z)-log(pdf(pd_cost_per_case_reduction_5,z)),(x0),[],[],[],[],(lb),(ub));
    Cost_Difference{3} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_cost_per_case_reduction_5,0.025),'%4.2f') char(8211) num2str(icdf(pd_cost_per_case_reduction_5,0.975),'%4.2f') ')'];
    Cost_Difference_RAW{3} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(prctile(DX_5(:),2.5),'%4.2f') char(8211) num2str(prctile(DX_5(:),97.5),'%4.2f') ')'];
Cost_per_Case_Difference=table(Scenario,Cost_Difference,Cost_Difference_RAW);
end