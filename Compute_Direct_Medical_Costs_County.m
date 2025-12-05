function [Cost_Case_Medical]= Compute_Direct_Medical_Costs_County(payer_type,County_Data_Vaccine_Reduction,PT_Unvaccinated_Cases_County,PT_Vaccinated_Cases_County,p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization,Cost_per_Non_Hospitalization)


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
        Adjustment_Cost_Inpatient=medc_ratio_in.*[mm_ratio.*ones(1,13) ones(1,5)];

        medc_ratio_out=Medicare_Ratio.Outpatient(f_medcare)./100;
        Adjustment_Cost_Outpatient=medc_ratio_out.*[mm_ratio.*ones(1,13) ones(1,5)];
    elseif strcmp(payer_type,'Private')   
        Adjustment_Cost_Inpatient=[ones(1,18)];
        Adjustment_Cost_Outpatient=[ones(1,18)];
    end

    Cost_Hospitalizations(ww,:,:)=State_Hospital_Cost.Cost_2025(f_state).*repmat(Adjustment_Cost_Inpatient(:),1,size(PT_Unvaccinated_Cases_County,3)).*(repmat(p_H_Unvaccinated(:).*duration_hospitalization(:),1,size(PT_Unvaccinated_Cases_County,3)).*squeeze(PT_Unvaccinated_Cases_County(ww,:,:))+repmat(p_H_Vaccinated(:).*duration_hospitalization(:),1,size(PT_Vaccinated_Cases_County,3)).*squeeze(PT_Vaccinated_Cases_County(ww,:,:)));    
    Cost_Non_Hospitalizations(ww,:,:)=Cost_per_Non_Hospitalization.*repmat(Adjustment_Cost_Outpatient(:),1,size(PT_Unvaccinated_Cases_County,3)).*(repmat(1-p_H_Unvaccinated(:),1,size(PT_Unvaccinated_Cases_County,3)).*squeeze(PT_Unvaccinated_Cases_County(ww,:,:))+repmat(1-p_H_Vaccinated(:),1,size(PT_Vaccinated_Cases_County,3)).*squeeze(PT_Vaccinated_Cases_County(ww,:,:)));

end


Cost_Case_Medical=squeeze(sum(Cost_Hospitalizations,[2]))+squeeze(sum(Cost_Non_Hospitalizations,[2]));
end