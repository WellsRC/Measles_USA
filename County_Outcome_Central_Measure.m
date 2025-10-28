function [cases,hospital,cost,cost_per_case]=County_Outcome_Central_Measure(National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6)

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
[p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();
[Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

[Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,2));

Cost_Vaccination_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
for cc=1:size(Hospitalizations_Baseline,1)
    Hospitalizations_Baseline(cc,:)=p_H_Unvaccinated*squeeze(Unvaccinated_Cases_County_Baseline(cc,:,:))+p_H_Vaccinated*squeeze(Vaccinated_Cases_County_Baseline(cc,:,:));
    Cost_Vaccination_Contacts(cc,:)=Cost_Vac_Age*squeeze(Unvaccinated_Contacts_Baseline(cc,:,:));
end

Cost_Case_Medical=Cost_per_Hospitalization.*Hospitalizations_Baseline+Cost_per_Non_Hospitalization.*(Outbreak_Cases_County-Hospitalizations_Baseline);
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;
Cost_Baseline=Contact_Tracing_Costs+Cost_Case_Medical+Testing_Cost+Cost_Vaccination_Contacts;

cases=mean(Outbreak_Cases_County,2);
hospital=mean(Hospitalizations_Baseline,2);

cost=NaN.*zeros(size(hospital));
cost_per_case=NaN.*zeros(size(hospital));
for cc=1:length(cost_per_case)
    tf=Outbreak_Cases_County(cc,:)>0;
    if(sum(tf)>0)
        cost(cc)=mean(Cost_Baseline(cc,tf));
        cost_per_case(cc)=mean(Cost_Baseline(cc,tf)./Outbreak_Cases_County(cc,tf));
    end
end

end
