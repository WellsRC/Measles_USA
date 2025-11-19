function [Total_Cost,Total_Outbreak_Response,Total_Medical,Total_Productivity_loss]=County_Level_Costs(County_Data,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline)

[p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization]=Hospitalization_Probability();
[Productivity_Days_Lost_Under_15_Case,Productivity_Days_Lost_15_plus_Case,Productivity_Days_Lost_Under_15_Contact,Productivity_Days_Lost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();


Productivity_loss_Cases=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) Productivity_Days_Lost_15_plus_Case.*ones(1,15)];
Productivity_loss_Cases_School=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) zeros(1,15)];

Productivity_loss_Contacts=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3)./3 Productivity_Days_Lost_15_plus_Contact.*ones(1,15)]; % Divide by 3 as assume productivity loss is 1 caregivr to 3 children
Productivity_loss_Contacts_School=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3) zeros(1,15)]; 

School_Cost=readtable('Cost_School_Day.xlsx');
State_Hospital_Cost=readtable('State_Daily_Hospital_Cost.xlsx');
Wages=readtable('County_Daily_Wages_2025.xlsx');
W=zeros(length(County_Data.GEOID),1);
Cost_Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County_Baseline));
Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));
Cost_Non_Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));

Total_Productivity_loss_Cases=zeros(size(Unvaccinated_Cases_County_Baseline));
Total_Productivity_loss_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline));


% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];

Cost_Vaccination_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline,1),size(Unvaccinated_Cases_County_Baseline,3));


for ww=1:length(W)
    f_state=strcmp(State_Hospital_Cost.State,County_Data.State{ww});
    school_state=strcmp(School_Cost.State,County_Data.State{ww});
    f_wages=strcmp(Wages.area_fips,County_Data.GEOID{ww});
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Cases(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Cases_School(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    Total_Productivity_loss_Cases(ww,:,:)=Wage_Loss+School_Loss;
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Contacts(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Contacts_School(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    Total_Productivity_loss_Contacts(ww,:,:)=Wage_Loss+School_Loss;
    Cost_Hospitalizations_Baseline(ww,:,:)=State_Hospital_Cost.Cost_2025(f_state).*(repmat(p_H_Unvaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:))+repmat(p_H_Vaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Vaccinated_Cases_County_Baseline(ww,:,:)));    
    Hospitalizations_Baseline(ww,:)=p_H_Unvaccinated*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:))+p_H_Vaccinated*squeeze(Vaccinated_Cases_County_Baseline(ww,:,:));
    Cost_Non_Hospitalizations_Baseline(ww,:,:)=Cost_per_Non_Hospitalization.*((1-p_H_Unvaccinated)*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)))+(1-p_H_Vaccinated)*squeeze(Vaccinated_Cases_County_Baseline(ww,:,:));
    Cost_Vaccination_Contacts(ww,:)=Cost_Vac_Age*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
end

Total_Medical=squeeze(sum(Cost_Hospitalizations_Baseline,[2]))+Cost_Non_Hospitalizations_Baseline;
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*squeeze(sum(Total_Contacts_Baseline,2));
Contact_Tracing_Costs=Cost_per_Contact.*squeeze(sum(Total_Contacts_Baseline,2));

Total_Productivity_loss_Cases=squeeze(sum(Total_Productivity_loss_Cases,2));
Total_Productivity_loss_Contacts=squeeze(sum(Total_Productivity_loss_Contacts,2));

Total_Productivity_loss=Total_Productivity_loss_Cases+Total_Productivity_loss_Contacts;
Total_Outbreak_Response=Contact_Tracing_Costs+Testing_Cost+Cost_Vaccination_Contacts;

Total_Cost=Contact_Tracing_Costs+Total_Medical+Testing_Cost+Cost_Vaccination_Contacts+Total_Productivity_loss_Cases+Total_Productivity_loss_Contacts;

end