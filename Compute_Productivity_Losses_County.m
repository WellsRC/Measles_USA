function [Total_Productivity_loss_Cases,Total_Productivity_loss_Contacts]=Compute_Productivity_Losses_County(County_Data_Vaccine_Reduction,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Unvaccinated_Contacts_Baseline,Productivity_Days_Lost_Under_15_Case,Productivity_Days_Lost_15_plus_Case,Productivity_Days_Lost_Under_15_Contact,Productivity_Days_Lost_15_plus_Contact)

Productivity_loss_Cases=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) Productivity_Days_Lost_15_plus_Case.*ones(1,15)];
Productivity_loss_Cases_School=[Productivity_Days_Lost_Under_15_Case.*ones(1,3) zeros(1,15)];

Productivity_loss_Contacts=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3)./3 Productivity_Days_Lost_15_plus_Contact.*ones(1,15)]; % Divide by 3 as assume productivity loss is 1 caregivr to 3 children
Productivity_loss_Contacts_School=[Productivity_Days_Lost_Under_15_Contact.*ones(1,3) zeros(1,15)]; 

School_Cost=readtable('Cost_School_Day.xlsx');

Wages=readtable('County_Daily_Wages_2025.xlsx');



Total_Productivity_loss_Cases=zeros(size(Unvaccinated_Cases_County_Baseline));
Total_Productivity_loss_Contacts=zeros(size(Unvaccinated_Cases_County_Baseline));


for ww=1:length(County_Data_Vaccine_Reduction.GEOID)
    school_state=strcmp(School_Cost.State,County_Data_Vaccine_Reduction.State{ww});
    f_wages=strcmp(Wages.area_fips,County_Data_Vaccine_Reduction.GEOID{ww});
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Cases(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Cases_School(:),1,size(Unvaccinated_Cases_County_Baseline,3)).*squeeze(Unvaccinated_Cases_County_Baseline(ww,:,:)+Vaccinated_Cases_County_Baseline(ww,:,:));
    Total_Productivity_loss_Cases(ww,:,:)=Wage_Loss+School_Loss;
    Wage_Loss=Wages.daily_wage_2025(f_wages).*repmat(Productivity_loss_Contacts(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    School_Loss=School_Cost.Daily_Cost_Student(school_state).*repmat(Productivity_loss_Contacts_School(:),1,size(Unvaccinated_Contacts_Baseline,3)).*squeeze(Unvaccinated_Contacts_Baseline(ww,:,:));
    Total_Productivity_loss_Contacts(ww,:,:)=Wage_Loss+School_Loss;           
end

Total_Productivity_loss_Cases=squeeze(sum(Total_Productivity_loss_Cases,[2]));
Total_Productivity_loss_Contacts=squeeze(sum(Total_Productivity_loss_Contacts,[2]));
end