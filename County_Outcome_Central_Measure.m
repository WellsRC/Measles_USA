function [cases,hospital,cost,cost_per_case]=County_Outcome_Central_Measure(National_Reduction,Scenario_Plot,Age_0_to_6)

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
[p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();
[Productivity_Cost_Under_15_Case,Productivity_Cost_15_plus_Case,Productivity_Cost_Under_15_Contact,Productivity_Cost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

if(Age_0_to_6)
    load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(National_Reduction.*100) '_Ages_0_to_6.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');
else
    load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(National_Reduction.*100) '_Ages_0_to_4.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');
end

Productivity_loss_Cases=[Productivity_Cost_Under_15_Case.*ones(1,3) Productivity_Cost_15_plus_Case.*ones(1,15)];
Productivity_loss_Contacts=[Productivity_Cost_Under_15_Contact.*ones(1,3) Productivity_Cost_15_plus_Contact.*ones(1,15)];

% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,2));

Cost_Vaccination_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
Productivity_Loss_total=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
for cc=1:size(Hospitalizations_Baseline,1)
    Hospitalizations_Baseline(cc,:)=p_H_Unvaccinated*squeeze(Unvaccinated_Cases_County_Baseline(cc,:,:))+p_H_Vaccinated*squeeze(Vaccinated_Cases_County_Baseline(cc,:,:));
    Cost_Vaccination_Contacts(cc,:)=Cost_Vac_Age*squeeze(Unvaccinated_Contacts_Baseline(cc,:,:));

    Productivity_Loss_total(cc,:)=Productivity_loss_Contacts*squeeze(Unvaccinated_Contacts_Baseline(cc,:,:))+Productivity_loss_Cases*(squeeze(Unvaccinated_Cases_County_Baseline(cc,:,:))+squeeze(Vaccinated_Cases_County_Baseline(cc,:,:)));
end

Cost_Case_Medical=Cost_per_Hospitalization.*Hospitalizations_Baseline+Cost_per_Non_Hospitalization.*(Total_Cases_County-Hospitalizations_Baseline);
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;
Cost_Baseline=Contact_Tracing_Costs+Cost_Case_Medical+Testing_Cost+Cost_Vaccination_Contacts+Productivity_Loss_total;

cases=mean(Total_Cases_County,2);
hospital=mean(Hospitalizations_Baseline,2);

cost=NaN.*zeros(size(hospital));
cost_per_case=NaN.*zeros(size(hospital));
for cc=1:length(cost_per_case)
    tf=Total_Cases_County(cc,:)>0;
    if(sum(tf)>0)
        cost(cc)=mean(Cost_Baseline(cc,tf));
        cost_per_case(cc)=mean(Cost_Baseline(cc,tf)./Total_Cases_County(cc,tf));
    end
end

end
