function [pd_cases,pd_hospital,pd_cost]=National_Outcome_Comparison_Distribution(National_Annual_Reduction,Scenario_Plot,Year_Reduced)

[p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization]=Hospitalization_Probability();
[Productivity_Days_Lost_Under_15_Case,Productivity_Days_Lost_15_plus_Case,Productivity_Days_Lost_Under_15_Contact,Productivity_Days_Lost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost();

    load(['National_Reduction=' num2str(100.*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'County_Data_Vaccine_Reduction')
    load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(100.*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');

Productivity_loss_Cases=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) Productivity_Days_Lost_15_plus_Case.*ones(1,15)];
Productivity_loss_Cases_School=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) zeros(1,15)];

Productivity_loss_Contacts=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3)./3 Productivity_Days_Lost_15_plus_Contact.*ones(1,15)]; % Divide by 3 as assume productivity loss is 1 caregivr to 3 children
Productivity_loss_Contacts_School=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3) zeros(1,15)]; 

School_Cost=readtable('Cost_School_Day.xlsx');
State_Hospital_Cost=readtable('State_Daily_Hospital_Cost.xlsx');
Wages=readtable('County_Daily_Wages_2025.xlsx');
W=zeros(length(County_Data_Vaccine_Reduction.GEOID),1);
Cost_Hospitalizations=zeros(size(Unvaccinated_Cases_County_Baseline));

Total_Productivity_loss_Cases=zeros(size(Unvaccinated_Cases_County_Baseline));
Total_Productivity_loss_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline));


Hospitalizations_Reduction=p_H_Unvaccinated*squeeze(sum(Unvaccinated_Cases_County_Baseline,1))+p_H_Vaccinated*squeeze(sum(Vaccinated_Cases_County_Baseline,1));
Cost_Non_Hospitalizations=Cost_per_Non_Hospitalization.*((1-p_H_Unvaccinated)*squeeze(sum(Unvaccinated_Cases_County_Baseline,1))+(1-p_H_Vaccinated)*squeeze(sum(Vaccinated_Cases_County_Baseline,1)));

for ww=1:length(W)
    f_state=strcmp(State_Hospital_Cost.State,County_Data_Vaccine_Reduction.State{ww});
    school_state=strcmp(School_Cost.State,County_Data_Vaccine_Reduction.State{ww});
    f_wages=strcmp(Wages.area_fips,County_Data_Vaccine_Reduction.GEOID{ww});
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Cases(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Cases_School(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    Total_Productivity_loss_Cases(ww,:,:)=Wage_Loss+School_Loss;
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Contacts(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Contacts_School(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    Total_Productivity_loss_Contacts(ww,:,:)=Wage_Loss+School_Loss;
    Cost_Hospitalizations(ww,:,:)=State_Hospital_Cost.Cost_2025(f_state).*(repmat(p_H_Unvaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:))+repmat(p_H_Vaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Vaccinated_Cases_County_Baseline(ww,:,:)));    
end



Total_Productivity_loss_Cases=sum(Total_Productivity_loss_Cases,[1 2]);
Total_Productivity_loss_Contacts=sum(Total_Productivity_loss_Contacts,[1 2]);

Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,[1 2]));
% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Cost_Vaccination_Contacts=Cost_Vac_Age*squeeze(sum(Unvaccinated_Contacts_Baseline,1));
Cost_Case_Medical=squeeze(sum(Cost_Hospitalizations,[1 2]))+Cost_Non_Hospitalizations(:);
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;

Contact_Tracing_Costs=Contact_Tracing_Costs(:);


Total_Cases_County_Reduction=Total_Cases_County;
Cost_Reduction=Contact_Tracing_Costs(:)+Cost_Case_Medical(:)+Testing_Cost(:)+Cost_Vaccination_Contacts(:)+Total_Productivity_loss_Cases(:)+Total_Productivity_loss_Contacts(:);

    load(['National_Reduction=' num2str(0) '_Year=' num2str(0) '.mat'],'County_Data_Vaccine_Reduction')
    load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(0) '_Year=' num2str(0) '.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');


Productivity_loss_Cases=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) Productivity_Days_Lost_15_plus_Case.*ones(1,15)];
Productivity_loss_Cases_School=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) zeros(1,15)];

Productivity_loss_Contacts=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3)./3 Productivity_Days_Lost_15_plus_Contact.*ones(1,15)]; % Divide by 3 as assume productivity loss is 1 caregivr to 3 children
Productivity_loss_Contacts_School=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3) zeros(1,15)]; 

School_Cost=readtable('Cost_School_Day.xlsx');
State_Hospital_Cost=readtable('State_Daily_Hospital_Cost.xlsx');
Wages=readtable('County_Daily_Wages_2025.xlsx');
W=zeros(length(County_Data_Vaccine_Reduction.GEOID),1);
Cost_Hospitalizations=zeros(size(Unvaccinated_Cases_County_Baseline));

Total_Productivity_loss_Cases=zeros(size(Unvaccinated_Cases_County_Baseline));
Total_Productivity_loss_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline));


Hospitalizations_Baseline=p_H_Unvaccinated*squeeze(sum(Unvaccinated_Cases_County_Baseline,1))+p_H_Vaccinated*squeeze(sum(Vaccinated_Cases_County_Baseline,1));
Cost_Non_Hospitalizations=Cost_per_Non_Hospitalization.*((1-p_H_Unvaccinated)*squeeze(sum(Unvaccinated_Cases_County_Baseline,1))+(1-p_H_Vaccinated)*squeeze(sum(Vaccinated_Cases_County_Baseline,1)));

for ww=1:length(W)
    f_state=strcmp(State_Hospital_Cost.State,County_Data_Vaccine_Reduction.State{ww});
    school_state=strcmp(School_Cost.State,County_Data_Vaccine_Reduction.State{ww});
    f_wages=strcmp(Wages.area_fips,County_Data_Vaccine_Reduction.GEOID{ww});
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Cases(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Cases_School(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    Total_Productivity_loss_Cases(ww,:,:)=Wage_Loss+School_Loss;
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Contacts(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Contacts_School(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    Total_Productivity_loss_Contacts(ww,:,:)=Wage_Loss+School_Loss;
    Cost_Hospitalizations(ww,:,:)=State_Hospital_Cost.Cost_2025(f_state).*(repmat(p_H_Unvaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:))+repmat(p_H_Vaccinated(:).*duration_hospitalization(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Vaccinated_Cases_County_Baseline(ww,:,:)));    
end



Total_Productivity_loss_Cases=sum(Total_Productivity_loss_Cases,[1 2]);
Total_Productivity_loss_Contacts=sum(Total_Productivity_loss_Contacts,[1 2]);

Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,[1 2]));
% https://pmc.ncbi.nlm.nih.gov/articles/PMC11309373/#:~:text=In%202023%2C%20approximately%2054%25%20of,)%20born%20during%201994%E2%80%932023.
% 54 % elgible for VFC
Cost_Vac_Age=[(0.54.*Cost_per_Vaccine_dose_VFC+(1-0.54).*Cost_per_Vaccine_dose_Private).*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
Cost_Vaccination_Contacts=Cost_Vac_Age*squeeze(sum(Unvaccinated_Contacts_Baseline,1));
Cost_Case_Medical=squeeze(sum(Cost_Hospitalizations,[1 2]))+Cost_Non_Hospitalizations(:);
Testing_Cost=Tests_per_Contact.*Cost_per_Test.*Total_Contacts_Baseline;
Contact_Tracing_Costs=Cost_per_Contact.*Total_Contacts_Baseline;

Contact_Tracing_Costs=Contact_Tracing_Costs(:);

Cost_Baseline=Contact_Tracing_Costs(:)+Cost_Case_Medical(:)+Testing_Cost(:)+Cost_Vaccination_Contacts(:)+Total_Productivity_loss_Cases(:)+Total_Productivity_loss_Contacts(:);



% Cases    
temp_c=sum(Total_Cases_County_Reduction-Total_Cases_County,1);
if(min(temp_c)>=0)
    temp_c(temp_c==0)=10^(-8);
    pd_cases=fitdist(temp_c(:),'Kernel','Support','positive');    
else
    pd_cases=fitdist(temp_c(:),'Kernel');
end
% Hospital
temp_c=Hospitalizations_Reduction-Hospitalizations_Baseline;
if(min(temp_c)>=0)
    temp_c(temp_c==0)=10^(-8);
    pd_hospital=fitdist(temp_c(:),'Kernel','Support','positive');    
else
    pd_hospital=fitdist(temp_c(:),'Kernel');
end

% Cost
temp_c=(Cost_Reduction-Cost_Baseline)./10^6;
if(min(temp_c)>=0)
    temp_c(temp_c==0)=10^(-8);
    pd_cost=fitdist(temp_c(:),'Kernel','Support','positive');
else
    pd_cost=fitdist(temp_c(:),'Kernel');
end

end
