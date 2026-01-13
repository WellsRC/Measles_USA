clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=1;

Output_Table=cell(4,1);
Cumulative_Costs=cell(5,4);
for National_Annual_Reduction=0.0025:0.0025:0.01
    [pd_1_cases,pd_1_hospital,pd_1_cost,pd_1_cost_per_case,pd_1_pro_loss,pd_1_med_cost,pd_1_med_cost_uninsured,pd_1_med_cost_public,pd_1_med_cost_private,~,~,pd_1_outbreak_response_cost,~,pd_1_deaths]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',1);
    [pd_2_cases,pd_2_hospital,pd_2_cost,pd_2_cost_per_case,pd_2_pro_loss,pd_2_med_cost,pd_2_med_cost_uninsured,pd_2_med_cost_public,pd_2_med_cost_private,~,~,pd_2_outbreak_response_cost,~,pd_2_deaths]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',2);
    [pd_3_cases,pd_3_hospital,pd_3_cost,pd_3_cost_per_case,pd_3_pro_loss,pd_3_med_cost,pd_3_med_cost_uninsured,pd_3_med_cost_public,pd_3_med_cost_private,~,~,pd_3_outbreak_response_cost,~,pd_3_deaths]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',3);
    [pd_4_cases,pd_4_hospital,pd_4_cost,pd_4_cost_per_case,pd_4_pro_loss,pd_4_med_cost,pd_4_med_cost_uninsured,pd_4_med_cost_public,pd_4_med_cost_private,~,~,pd_4_outbreak_response_cost,~,pd_4_deaths]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',4);
    [pd_5_cases,pd_5_hospital,pd_5_cost,pd_5_cost_per_case,pd_5_pro_loss,pd_5_med_cost,pd_5_med_cost_uninsured,pd_5_med_cost_public,pd_5_med_cost_private,~,~,pd_5_outbreak_response_cost,~,pd_5_deaths]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',5);
     
    Year={'Y1';'Y2';'Y3';'Y4';'Y5'};
    Cases=cell(length(Year),1);
    Hospitalizations=cell(length(Year),1);
    Deaths=cell(length(Year),1);
    Cost=cell(length(Year),1);
    Cost_per_Case=cell(length(Year),1);
    Outbreak_Response=cell(length(Year),1);
    Direct_Medical=cell(length(Year),1);
    Direct_Medical_Uninsured=cell(length(Year),1);
    Direct_Medical_Public=cell(length(Year),1);
    Direct_Medical_Private=cell(length(Year),1);
    Productivity_Loss=cell(length(Year),1);

    % Cases
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_cases,0.5);    
    Cases{1} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_cases,0.5); 
    Cases{2} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_cases,0.5); 
    Cases{3} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_cases,0.5); 
    Cases{4} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_cases,0.5); 
    Cases{5} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    % Hospitalizations
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_hospital,0.5);    
    Hospitalizations{1} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_hospital,0.5); 
    Hospitalizations{2} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_hospital,0.5); 
    Hospitalizations{3} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_hospital,0.5); 
    Hospitalizations{4} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_hospital,0.5); 
    Hospitalizations{5} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];

    % Deaths
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_deaths,0.5);    
    Deaths{1} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_deaths,0.5); 
    Deaths{2} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_deaths,0.5); 
    Deaths{3} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_deaths,0.5); 
    Deaths{4} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_deaths,0.5); 
    Deaths{5} = [num2str(mle,'%5.0f') ' (' num2str(lb_hdi,'%5.0f') char(8211) num2str(ub_hdi,'%5.0f') ')'];
    
    % Cost
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_cost,0.5);    
    Cost{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_cost,0.5); 
    Cost{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_cost,0.5); 
    Cost{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_cost,0.5); 
    Cost{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_cost,0.5); 
    Cost{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    
    % Outbreak_repsponse

    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_outbreak_response_cost,0.5);    
    Outbreak_Response{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_outbreak_response_cost,0.5); 
    Outbreak_Response{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_outbreak_response_cost,0.5); 
    Outbreak_Response{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_outbreak_response_cost,0.5); 
    Outbreak_Response{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_outbreak_response_cost,0.5); 
    Outbreak_Response{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    
    % Direct Medical
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_med_cost,0.5);    
    Direct_Medical{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_med_cost,0.5); 
    Direct_Medical{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_med_cost,0.5); 
    Direct_Medical{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_med_cost,0.5); 
    Direct_Medical{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_med_cost,0.5); 
    Direct_Medical{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    

     % Direct Medical: Uninsured

     [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_med_cost_uninsured,0.5);    
    Direct_Medical_Uninsured{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_med_cost_uninsured,0.5); 
    Direct_Medical_Uninsured{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_med_cost_uninsured,0.5); 
    Direct_Medical_Uninsured{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_med_cost_uninsured,0.5); 
    Direct_Medical_Uninsured{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_med_cost_uninsured,0.5); 
    Direct_Medical_Uninsured{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    % Direct Medical: Public

    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_med_cost_public,0.5);    
    Direct_Medical_Public{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_med_cost_public,0.5); 
    Direct_Medical_Public{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_med_cost_public,0.5); 
    Direct_Medical_Public{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_med_cost_public,0.5); 
    Direct_Medical_Public{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_med_cost_public,0.5); 
    Direct_Medical_Public{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    % Direct Medical: Private
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_med_cost_private,0.5);    
    Direct_Medical_Private{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_med_cost_private,0.5); 
    Direct_Medical_Private{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_med_cost_private,0.5); 
    Direct_Medical_Private{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_med_cost_private,0.5); 
    Direct_Medical_Private{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_med_cost_private,0.5); 
    Direct_Medical_Private{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];

    % Productivity loss

    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_pro_loss,0.5);    
    Productivity_Loss{1} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_pro_loss,0.5); 
    Productivity_Loss{2} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_pro_loss,0.5); 
    Productivity_Loss{3} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_pro_loss,0.5); 
    Productivity_Loss{4} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_pro_loss,0.5); 
    Productivity_Loss{5} = [num2str(mle,'%4.1f') ' (' num2str(lb_hdi,'%4.1f') char(8211) num2str(ub_hdi,'%4.1f') ')'];

    % Cost per case

    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_cost_per_case,0.5);    
    Cost_per_Case{1} = [num2str(mle,'%6.0f') ' (' num2str(lb_hdi,'%6.0f') char(8211) num2str(ub_hdi,'%6.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_2_cost_per_case,0.5); 
    Cost_per_Case{2} = [num2str(mle,'%6.0f') ' (' num2str(lb_hdi,'%6.0f') char(8211) num2str(ub_hdi,'%6.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_3_cost_per_case,0.5); 
    Cost_per_Case{3} = [num2str(mle,'%6.0f') ' (' num2str(lb_hdi,'%6.0f') char(8211) num2str(ub_hdi,'%6.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_4_cost_per_case,0.5); 
    Cost_per_Case{4} = [num2str(mle,'%6.0f') ' (' num2str(lb_hdi,'%6.0f') char(8211) num2str(ub_hdi,'%6.0f') ')'];
    
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_5_cost_per_case,0.5); 
    Cost_per_Case{5} = [num2str(mle,'%6.0f') ' (' num2str(lb_hdi,'%6.0f') char(8211) num2str(ub_hdi,'%6.0f') ')'];
    
    Output_Table{count}=table(Year,Cases,Hospitalizations,Deaths, Cost,Outbreak_Response,Direct_Medical,Direct_Medical_Uninsured,Direct_Medical_Public,Direct_Medical_Private,Productivity_Loss,Cost_per_Case);
    

    rng(30302030);
    r1=random(pd_1_cost,2500,1);
    r2=random(pd_2_cost,2500,1);
    r3=random(pd_3_cost,2500,1);
    r4=random(pd_4_cost,2500,1);
    r5=random(pd_5_cost,2500,1);

    
    % T=r1;
    % pdC=fitdist(T./10^3,'Kernel','Support','positive');
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pd_1_cost,0.5); 
    Cumulative_Costs{1,count}=['$' num2str(mle./10^3,'%3.2f') ' ($' num2str(lb_hdi./10^3,'%3.2f') char(8211) '$' num2str(ub_hdi./10^3,'%3.2f') ')'];

    T=r1+r2;
    pdC=fitdist(T./10^3,'Kernel','Support','positive');
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pdC,0.5); 
    Cumulative_Costs{2,count}=['$' num2str(mle,'%3.2f') ' ($' num2str(lb_hdi,'%3.2f') char(8211) '$' num2str(ub_hdi,'%3.2f') ')'];
    
    T=r1+r2+r3;
    pdC=fitdist(T./10^3,'Kernel','Support','positive');
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pdC,0.5);
    Cumulative_Costs{3,count}=['$' num2str(mle,'%3.2f') ' ($' num2str(lb_hdi,'%3.2f') char(8211) '$' num2str(ub_hdi,'%3.2f') ')'];

    T=r1+r2+r3+r4;
    pdC=fitdist(T./10^3,'Kernel','Support','positive');
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pdC,0.5);
    Cumulative_Costs{4,count}=['$' num2str(mle,'%3.2f') ' ($' num2str(lb_hdi,'%3.2f') char(8211) '$' num2str(ub_hdi,'%3.2f') ')'];
    
    T=r1+r2+r3+r4+r5;
    pdC=fitdist(T./10^3,'Kernel','Support','positive');
    [mle,lb_hdi,ub_hdi]=Estimate_HDI(pdC,0.5);
    Cumulative_Costs{5,count}=['$' num2str(mle,'%3.2f') ' ($' num2str(lb_hdi,'%3.2f') char(8211) '$' num2str(ub_hdi,'%3.2f') ')'];
    
    count=count+1;
end