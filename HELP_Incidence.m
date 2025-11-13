clear;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[pd_cases,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_test_vac_cost,pd_ct_cost,pd_pro_loss_per_case,pd_med_cost_per_case,pd_test_vac_cost_per_case,pd_ct_cost_per_case,break_pro_loss,break_test_vac_cost,break_ct_cost,break_med_cost]=National_Outcome_Distribution(0,'Sample_2025',0);

National_Annual_Reduction=0.01;
[pd_1_cases,pd_1_hospital,pd_1_cost]=National_Outcome_Distribution(National_Annual_Reduction,'Sample_2025',1);
[pd_2_cases,pd_2_hospital,pd_2_cost]=National_Outcome_Distribution(National_Annual_Reduction,'Sample_2025',2);
[pd_3_cases,pd_3_hospital,pd_3_cost]=National_Outcome_Distribution(National_Annual_Reduction,'Sample_2025',3);
[pd_4_cases,pd_4_hospital,pd_4_cost]=National_Outcome_Distribution(National_Annual_Reduction,'Sample_2025',4);
[pd_5_cases,pd_5_hospital,pd_5_cost]=National_Outcome_Distribution(National_Annual_Reduction,'Sample_2025',5);
close all;

Year={'Baseline';'1% reduction';'2% reduction';'3% reduction';'4% reduction';'5% reduction'};
Cases=cell(length(Year),1);
Hospitalizations=cell(length(Year),1);
Cost=cell(length(Year),1);


% Cases
x0=icdf(pd_cases,0.5);
lb=icdf(pd_cases,0.05);
ub=icdf(pd_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{1} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_cases,0.975),'%5.0f') ')'];

x0=icdf(pd_1_cases,0.5);
lb=icdf(pd_1_cases,0.05);
ub=icdf(pd_1_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_1_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{2} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_1_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_1_cases,0.975),'%5.0f') ')'];

x0=icdf(pd_2_cases,0.5);
lb=icdf(pd_2_cases,0.05);
ub=icdf(pd_2_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_2_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{3} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_2_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_2_cases,0.975),'%5.0f') ')'];

x0=icdf(pd_3_cases,0.5);
lb=icdf(pd_3_cases,0.05);
ub=icdf(pd_3_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_3_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{4} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_3_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_3_cases,0.975),'%5.0f') ')'];

x0=icdf(pd_4_cases,0.5);
lb=icdf(pd_4_cases,0.05);
ub=icdf(pd_4_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_4_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{5} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_4_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_4_cases,0.975),'%5.0f') ')'];

x0=icdf(pd_5_cases,0.5);
lb=icdf(pd_5_cases,0.05);
ub=icdf(pd_5_cases,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_5_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cases{6} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_5_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_5_cases,0.975),'%5.0f') ')'];

% Hospitalizations
x0=icdf(pd_hospital,0.5);
lb=icdf(pd_hospital,0.05);
ub=icdf(pd_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{1} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.975),'%5.0f') ')'];

x0=icdf(pd_1_hospital,0.5);
lb=icdf(pd_1_hospital,0.05);
ub=icdf(pd_1_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_1_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{2} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_1_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_1_hospital,0.975),'%5.0f') ')'];

x0=icdf(pd_2_hospital,0.5);
lb=icdf(pd_2_hospital,0.05);
ub=icdf(pd_2_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_2_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{3} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_2_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_2_hospital,0.975),'%5.0f') ')'];

x0=icdf(pd_3_hospital,0.5);
lb=icdf(pd_3_hospital,0.05);
ub=icdf(pd_3_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_3_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{4} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_3_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_3_hospital,0.975),'%5.0f') ')'];

x0=icdf(pd_4_hospital,0.5);
lb=icdf(pd_4_hospital,0.05);
ub=icdf(pd_4_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_4_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{5} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_4_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_4_hospital,0.975),'%5.0f') ')'];

x0=icdf(pd_5_hospital,0.5);
lb=icdf(pd_5_hospital,0.05);
ub=icdf(pd_5_hospital,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_5_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Hospitalizations{6} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_5_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_5_hospital,0.975),'%5.0f') ')'];

% Cost
x0=icdf(pd_cost,0.5);
lb=icdf(pd_cost,0.05);
ub=icdf(pd_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{1} = ['$' num2str(temp_X,'%5.1f') '(95% CrI: $' num2str(icdf(pd_cost,0.025),'%5.1f') char(8211) num2str(icdf(pd_cost,0.975),'%5.1f') ')'];

x0=icdf(pd_1_cost,0.5);
lb=icdf(pd_1_cost,0.05);
ub=icdf(pd_1_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_1_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{2} = ['$' num2str(temp_X,'%5.2f') '(95% CrI: $' num2str(icdf(pd_1_cost,0.025),'%5.2f') char(8211) '$' num2str(icdf(pd_1_cost,0.975),'%5.2f') ')'];

x0=icdf(pd_2_cost,0.5);
lb=icdf(pd_2_cost,0.05);
ub=icdf(pd_2_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_2_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{3} = ['$' num2str(temp_X,'%5.2f') '(95% CrI: $' num2str(icdf(pd_2_cost,0.025),'%5.2f') char(8211) '$' num2str(icdf(pd_2_cost,0.975),'%5.2f') ')'];

x0=icdf(pd_3_cost,0.5);
lb=icdf(pd_3_cost,0.05);
ub=icdf(pd_3_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_3_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{4} = ['$' num2str(temp_X,'%5.2f') '(95% CrI: $' num2str(icdf(pd_3_cost,0.025),'%5.2f') char(8211) '$' num2str(icdf(pd_3_cost,0.975),'%5.2f') ')'];

x0=icdf(pd_4_cost,0.5);
lb=icdf(pd_4_cost,0.05);
ub=icdf(pd_4_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_4_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{5} = ['$' num2str(temp_X,'%5.2f') '(95% CrI: $' num2str(icdf(pd_4_cost,0.025),'%5.2f') char(8211) '$' num2str(icdf(pd_4_cost,0.975),'%5.2f') ')'];

x0=icdf(pd_5_cost,0.5);
lb=icdf(pd_5_cost,0.05);
ub=icdf(pd_5_cost,0.95);
temp_X=10.^fmincon(@(z)-log(pdf(pd_5_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
Cost{6} = ['$' num2str(temp_X,'%5.2f') '(95% CrI: $' num2str(icdf(pd_5_cost,0.025),'%5.2f') char(8211) '$' num2str(icdf(pd_5_cost,0.975),'%5.2f') ')'];

Output_Table=table(Year,Cases,Hospitalizations,Cost);

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.075 0.15 0.87 0.8]);
cases=linspace(0,50000,10001);
plot(cases,pdf(pd_cases,cases),'color','k','LineWidth',2); hold on;
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cases(1):500:5000],'Fontsize',16)
xlabel(['Annual measles cases'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([500,cases(end)])
print(gcf,['HELP_Incidence_National.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.075 0.15 0.87 0.8]);
hospital=linspace(0,8000,10001);
plot(hospital,pdf(pd_hospital,hospital),'color','k','LineWidth',2); hold on;
set(gca,'LineWidth',2,'Tickdir','out','XTick',[hospital(1):50:hospital(end)],'Fontsize',16)
xlabel(['Annual measles hospitalizations'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([hospital(1),hospital(end)])
print(gcf,['HELP_Hospital_National.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.095 0.15 0.87 0.8]);
cost=linspace(0,3000,10001);
plot(cost,pdf(pd_cost,cost),'color','k','LineWidth',2); hold on;
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cost(1):20:cost(end)],'Fontsize',16)
xlabel(['Annual measles cost (Millions)'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([cost(1),cost(end)])
xtickformat('$%,.0f')
print(gcf,['HELP_Cost_National.png'],'-dpng','-r300');



figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.095 0.15 0.87 0.8]);
cost_per_case=linspace(80,120,10001);
plot(cost_per_case,pdf(pd_cost_per_case,cost_per_case),'color','k','LineWidth',2); hold on;
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cost_per_case(1):5:cost_per_case(end)],'Fontsize',16)
xlabel(['Cost per measles case (Thousands)'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([cost_per_case(1),cost_per_case(end)])
xtickformat('$%,.0f')
print(gcf,['HELP_Cost_per_Case_National.png'],'-dpng','-r300');

