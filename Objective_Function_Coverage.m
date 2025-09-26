function J = Objective_Function_Coverage(x,County_Data_model,State_Data,filter_county,filter_state)

[beta_x,beta_insurance,hyp_par,hyp_par_gam_a] = Parameters(x);
[Estimated_Vaccination_Coverage] = Vaccination_Coverage([County_Data_model.X County_Data_model.XI County_Data_model.X2],beta_x,beta_insurance,County_Data_model);

v_county=Estimated_Vaccination_Coverage.Overall;
v_state = zeros(height(State_Data),1);
v_state_Uninsured = zeros(height(State_Data),1);
v_state_Public = zeros(height(State_Data),1);
v_state_Private = zeros(height(State_Data),1);

u_year=unique(State_Data.Year);

for yy=1:length(u_year)
    tf=State_Data.Year==u_year(yy);
    year_County=County_Data_model.Year==u_year(yy);
    v_state(tf) = Estimated_State_Vaccine_Uptake(v_county(year_County),County_Data_model.Age_5_to_9(year_County),County_Data_model.State_FIP(year_County),State_Data.State_FIP(tf));
    v_state_Uninsured(tf) = Estimated_State_Vaccine_Uptake(Estimated_Vaccination_Coverage.Uninsured(year_County),County_Data_model.Uninsured(year_County).*County_Data_model.Age_5_to_9(year_County),County_Data_model.State_FIP(year_County),State_Data.State_FIP(tf));
    v_state_Public(tf) = Estimated_State_Vaccine_Uptake(Estimated_Vaccination_Coverage.Public(year_County),County_Data_model.Public(year_County).*County_Data_model.Age_5_to_9(year_County),County_Data_model.State_FIP(year_County),State_Data.State_FIP(tf));
    v_state_Private(tf) = Estimated_State_Vaccine_Uptake(Estimated_Vaccination_Coverage.Private(year_County),County_Data_model.Private(year_County).*County_Data_model.Age_5_to_9(year_County),County_Data_model.State_FIP(year_County),State_Data.State_FIP(tf));

end

% Compute the county-level coverage for all time and space
z_data=County_Data_model.Vaccine_Uptake(:);
z_data=log(z_data./(1-z_data));

z_model=v_county(:);
z_model=log(z_model./(1-z_model));

L_county=normpdf(z_data,z_model,hyp_par(1));
% L_county(County_Data.Vaccine_Uptake(:)==1)=v_county(County_Data.Vaccine_Uptake(:)==1);
% L_county(County_Data.Vaccine_Uptake(:)==0)=1-v_county(County_Data.Vaccine_Uptake(:)==0);

L_county=log(L_county(filter_county & ~isnan(County_Data_model.Vaccine_Uptake(:))));

L_state=zeros(height(State_Data),1);

OR_uninsured_private=v_state_Uninsured.*(1-v_state_Private)./(v_state_Private.*(1-v_state_Uninsured));
OR__medicare_private=v_state_Public.*(1-v_state_Private)./(v_state_Private.*(1-v_state_Public));

for ss=1:length(L_state)
    lbs=betainv(0.0001,State_Data.a_Beta_Parameters_Vaccine_Uptake(ss),State_Data.b_Beta_Parameters_Vaccine_Uptake(ss));
    ubs=betainv(0.9999,State_Data.a_Beta_Parameters_Vaccine_Uptake(ss),State_Data.b_Beta_Parameters_Vaccine_Uptake(ss));

    vx=linspace(lbs,ubs,1001);
    vx=vx(vx>0 & vx<1);
    dx=vx(2)-vx(1);
    test2=dx.*sum(betapdf(vx(2:end),State_Data.a_Beta_Parameters_Vaccine_Uptake(ss),State_Data.b_Beta_Parameters_Vaccine_Uptake(ss)).*normpdf(log(vx(2:end)./(1-vx(2:end))),log(v_state(ss)./(1-v_state(ss))),hyp_par(2)));


    L_state(ss)=log(test2);
    
end

hyp_par_gam_b=OR_uninsured_private./hyp_par_gam_a(1);
L_state_odd_uninsured_private=log(gampdf(State_Data.OR_Uninsured_Private,hyp_par_gam_a(1),hyp_par_gam_b));
L_state_odd_uninsured_private(isnan(L_state_odd_uninsured_private))=0;

L_state_odd_uninsured_private=L_state_odd_uninsured_private(filter_state);

hyp_par_gam_b=OR__medicare_private./hyp_par_gam_a(2);
L_state_odd_medicare_private=log(gampdf(State_Data.OR_Public_Private,hyp_par_gam_a(2),hyp_par_gam_b));
L_state_odd_medicare_private(isnan(L_state_odd_medicare_private))=0;
L_state_odd_medicare_private=L_state_odd_medicare_private(filter_state);



L_state=L_state(filter_state);

J=-sum(L_county(:)) -sum(L_state(:)) -sum(L_state_odd_medicare_private(:)) -sum(L_state_odd_uninsured_private(:));
end

