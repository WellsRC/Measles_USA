clear;
clc;

Measure={'Cases';'Hospitalizations';'Cost';'Outbreak respons';'Direct Medical';'Productivity loses';'Cost per case'};

[pd_cases,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_test_vac_cost,pd_ct_cost,pd_pro_loss_per_case,pd_med_cost_per_case,pd_test_vac_cost_per_case,pd_ct_cost_per_case,break_pro_loss,break_test_vac_cost,break_ct_cost,break_med_cost,pd_outbreak_response_cost]=National_Outcome_Distribution(0,'Baseline',0);

Value=cell(7,1);
  Value{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];
    
    Value{2} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];
    
    Value{3} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ') million'];
    
    Value{4} = [num2str(icdf(pd_outbreak_response_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_outbreak_response_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_outbreak_response_cost,0.75)./10^6,'%4.1f') ') million'];
    
    Value{5} = [num2str(icdf(pd_med_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_med_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_med_cost,0.75)./10^6,'%4.1f') ') million'];

    Value{6} = [num2str(icdf(pd_pro_loss,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_pro_loss,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_pro_loss,0.75)./10^6,'%4.1f') ') million'];

    Value{7} = [num2str(10^3.*icdf(pd_cost_per_case,0.5),'%5.0f') ' (IQR:' num2str(10^3.*icdf(pd_cost_per_case,0.25),'%5.0f') char(8211) num2str(10^3.*icdf(pd_cost_per_case,0.75),'%5.0f') ')'];

 Baseline_Output_Table=table(Measure,Value)




Measure={'Cases';'Hospitalizations';'Cost';'Outbreak respons';'Direct Medical';'Productivity loses';'Cost per case'};

[pd_cases,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_test_vac_cost,pd_ct_cost,pd_pro_loss_per_case,pd_med_cost_per_case,pd_test_vac_cost_per_case,pd_ct_cost_per_case,break_pro_loss,break_test_vac_cost,break_ct_cost,break_med_cost,pd_outbreak_response_cost]=National_Outcome_Distribution(0.01,'Baseline',5);

Value=cell(7,1);
  Value{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];
    
    Value{2} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];
    
    Value{3} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ') million'];
    
    Value{4} = [num2str(icdf(pd_outbreak_response_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_outbreak_response_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_outbreak_response_cost,0.75)./10^6,'%4.1f') ') million'];
    
    Value{5} = [num2str(icdf(pd_med_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_med_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_med_cost,0.75)./10^6,'%4.1f') ') million'];

    Value{6} = [num2str(icdf(pd_pro_loss,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_pro_loss,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_pro_loss,0.75)./10^6,'%4.1f') ') million'];

    Value{7} = [num2str(10^3.*icdf(pd_cost_per_case,0.5),'%5.0f') ' (IQR:' num2str(10^3.*icdf(pd_cost_per_case,0.25),'%5.0f') char(8211) num2str(10^3.*icdf(pd_cost_per_case,0.75),'%5.0f') ')'];

 Reduction_Output_Table=table(Measure,Value)