function [Cost_Case_Medical]= Compute_Direct_Medical_Costs(payer_type,County_Data_Vaccine_Reduction,Scenario_Plot,National_Annual_Reduction,Year_Reduced,p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization,Cost_per_Non_Hospitalization)


load(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(100*National_Annual_Reduction) '_Year=' num2str(Year_Reduced) '.mat'],'Uninsured_Unvaccinated_Cases_County','Uninsured_Vaccinated_Cases_County','Public_Unvaccinated_Cases_County','Public_Vaccinated_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline');

if strcmp(payer_type,'Uninsured')
    PT_Unvaccinated_Cases_County=Uninsured_Unvaccinated_Cases_County;
    PT_Vaccinated_Cases_County=Uninsured_Vaccinated_Cases_County;
elseif strcmp(payer_type,'Private')
    PT_Unvaccinated_Cases_County=Unvaccinated_Cases_County_Baseline-Uninsured_Unvaccinated_Cases_County-Public_Unvaccinated_Cases_County;
    PT_Vaccinated_Cases_County=Vaccinated_Cases_County_Baseline-Uninsured_Vaccinated_Cases_County-Public_Vaccinated_Cases_County;
elseif strcmp(payer_type,'Public')
    PT_Unvaccinated_Cases_County=Public_Unvaccinated_Cases_County;
    PT_Vaccinated_Cases_County=Public_Vaccinated_Cases_County;
end
clear Unvaccinated_Cases_County_Baseline Vaccinated_Cases_County_Baseline 'Uninsured_Unvaccinated_Cases_County' 'Uninsured_Vaccinated_Cases_County' 'Public_Unvaccinated_Cases_County' 'Public_Vaccinated_Cases_County'

State_Hospital_Cost=readtable('State_Daily_Hospital_Cost.xlsx');

% Medicare to Medicade Ratio
% https://www.healthaffairs.org/doi/full/10.1377/hlthaff.2020.00611
Medicaid_Medicare_Ratio=readtable('Medicaid-to-Medicare fee index by service type.xlsx');
% Medicare to private
Medicare_Ratio=readtable('Relative_Medicare_Prices.xlsx');


Cost_Hospitalizations=zeros(size(PT_Unvaccinated_Cases_County));

Cost_Non_Hospitalizations=zeros(size(PT_Unvaccinated_Cases_County));

for ww=1:length(County_Data_Vaccine_Reduction.GEOID)

    f_state=strcmp(State_Hospital_Cost.State,County_Data_Vaccine_Reduction.State{ww});
   
    if strcmp(payer_type,'Uninsured')
    % Uninsured relative to private
        Adjustment_Cost_Inpatient=1./[1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107];
        Adjustment_Cost_Outpatient=1./[1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107 1.107];
    elseif strcmp(payer_type,'Public')    
         f_medcare=strcmp(Medicare_Ratio.State,County_Data_Vaccine_Reduction.State{ww});
        f_medaid=strcmp(Medicaid_Medicare_Ratio.State,County_Data_Vaccine_Reduction.State{ww});
        medc_ratio_in=Medicare_Ratio.Inpatient(f_medcare)./100;
        mm_ratio=Medicaid_Medicare_Ratio.PrimaryCare(f_medaid);
        Adjustment_Cost_Inpatient=[mm_ratio.*ones(1,13) ones(1,5)]./medc_ratio_in;

        medc_ratio_out=Medicare_Ratio.Outpatient(f_medcare)./100;
        Adjustment_Cost_Outpatient=[mm_ratio.*ones(1,13) ones(1,5)]./medc_ratio_out;
    elseif strcmp(payer_type,'Private')   
        Adjustment_Cost_Inpatient=[ones(1,18)];
        Adjustment_Cost_Outpatient=[ones(1,18)];
    end

try
    Cost_Hospitalizations(ww,:,:)=State_Hospital_Cost.Cost_2025(f_state).*repmat(Adjustment_Cost_Inpatient(:),1,size(PT_Unvaccinated_Cases_County,3)).*(repmat(p_H_Unvaccinated(:).*duration_hospitalization(:),1,size(PT_Unvaccinated_Cases_County,3)).*squeeze(PT_Unvaccinated_Cases_County(ww,:,:))+repmat(p_H_Vaccinated(:).*duration_hospitalization(:),1,size(PT_Vaccinated_Cases_County,3)).*squeeze(PT_Vaccinated_Cases_County(ww,:,:)));    
    Cost_Non_Hospitalizations(ww,:,:)=Cost_per_Non_Hospitalization.*repmat(Adjustment_Cost_Outpatient(:),1,size(PT_Unvaccinated_Cases_County,3)).*(repmat(1-p_H_Unvaccinated(:),1,size(PT_Unvaccinated_Cases_County,3)).*squeeze(PT_Unvaccinated_Cases_County(ww,:,:))+repmat(1-p_H_Vaccinated(:),1,size(PT_Vaccinated_Cases_County,3)).*squeeze(PT_Vaccinated_Cases_County(ww,:,:)));
catch
    test=0;
end

end


Cost_Case_Medical=squeeze(sum(Cost_Hospitalizations,[1 2]))+squeeze(sum(Cost_Non_Hospitalizations,[1 2]));
end