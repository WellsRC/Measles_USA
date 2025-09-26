function J = Objective_Adjustment_Coverage(x,beta_x,beta_insurance,County_Data_model,known_county,state_fill,u_state_fill,State_Data)

dZ_County = zeros(size(known_county));

dZ_County(known_county)=x(1:sum(known_county));

for ss=1:length(u_state_fill)
    t_f=state_fill & County_Data_model.State_FIP==u_state_fill(ss);
    dZ_County(t_f)=x(sum(known_county)+ss);
end

[Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model,dZ_County);

v_county=Estimated_Vaccination_Coverage.Overall;
v_state = Estimated_State_Vaccine_Uptake(v_county,County_Data_model.Weight,County_Data_model.State_FIP,State_Data.State_FIP);
   
VAC_STAT=State_Data.Vaccine_Uptake;

VAC_STAT=VAC_STAT(ismember(State_Data.State_FIP,u_state_fill));
v_state=v_state(ismember(State_Data.State_FIP,u_state_fill));


Z_State=10^4.*(VAC_STAT-v_state);
Z_county=10^4.*(v_county(known_county)-County_Data_model.Vaccine_Uptake(known_county));
J = [Z_State;Z_county];
end