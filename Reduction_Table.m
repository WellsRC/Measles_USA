clear;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[pd_cases]=National_Outcome_Distribution(0,'Baseline',0);
[pd_cases_2024]=National_Outcome_Distribution(0,'Sample_2024',0);
[pd_cases_2023]=National_Outcome_Distribution(0,'Sample_2023',0);
% 
National_Annual_Reduction=0.01;
[pd_1_cases,pd_1_hospital,pd_1_cost,pd_1_cost_per_case]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',1);
[pd_2_cases,pd_2_hospital,pd_2_cost,pd_2_cost_per_case]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',2);
[pd_3_cases,pd_3_hospital,pd_3_cost,pd_3_cost_per_case]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',3);
[pd_4_cases,pd_4_hospital,pd_4_cost,pd_4_cost_per_case]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',4);
[pd_5_cases,pd_5_hospital,pd_5_cost,pd_5_cost_per_case]=National_Outcome_Distribution(National_Annual_Reduction,'Baseline',5);
close all;

Year={'Baseline';'1% reduction';'2% reduction';'3% reduction';'4% reduction';'5% reduction'};
Cases=cell(length(Year),1);
Hospitalizations=cell(length(Year),1);
Cost=cell(length(Year),1);
Cost_per_Case=cell(length(Year),1);


% Cases
Cases{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];

Cases{2} = [num2str(icdf(pd_1_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_1_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_1_cases,0.75),'%5.0f') ')'];

Cases{3} = [num2str(icdf(pd_2_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_2_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_2_cases,0.75),'%5.0f') ')'];

Cases{4} = [num2str(icdf(pd_3_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_3_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_3_cases,0.75),'%5.0f') ')'];

Cases{5} = [num2str(icdf(pd_4_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_4_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_4_cases,0.75),'%5.0f') ')'];

Cases{6} = [num2str(icdf(pd_5_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_5_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_5_cases,0.75),'%5.0f') ')'];

% Hospitalizations
Hospitalizations{1} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];

Hospitalizations{2} = [num2str(icdf(pd_1_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_1_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_1_hospital,0.75),'%5.0f') ')'];

Hospitalizations{3} = [num2str(icdf(pd_2_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_2_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_2_hospital,0.75),'%5.0f') ')'];

Hospitalizations{4} = [num2str(icdf(pd_3_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_3_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_3_hospital,0.75),'%5.0f') ')'];

Hospitalizations{5} = [num2str(icdf(pd_4_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_4_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_4_hospital,0.75),'%5.0f') ')'];

Hospitalizations{6} = [num2str(icdf(pd_5_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_5_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_5_hospital,0.75),'%5.0f') ')'];

% Cost
Cost{1} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ')'];

Cost{2} = [num2str(icdf(pd_1_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_1_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_1_cost,0.75),'%4.1f') ')'];

Cost{3} = [num2str(icdf(pd_2_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_2_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_2_cost,0.75),'%4.1f') ')'];

Cost{4} = [num2str(icdf(pd_3_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_3_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_3_cost,0.75),'%4.1f') ')'];

Cost{5} = [num2str(icdf(pd_4_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_4_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_4_cost,0.75),'%4.1f') ')'];

Cost{6} = [num2str(icdf(pd_5_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_5_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_5_cost,0.75),'%4.1f') ')'];


% Cost per case
Cost_per_Case{1} = [num2str(icdf(pd_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost_per_case,0.75),'%4.1f') ')'];

Cost_per_Case{2} = [num2str(icdf(pd_1_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_1_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_1_cost_per_case,0.75),'%4.1f') ')'];

Cost_per_Case{3} = [num2str(icdf(pd_2_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_2_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_2_cost_per_case,0.75),'%4.1f') ')'];

Cost_per_Case{4} = [num2str(icdf(pd_3_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_3_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_3_cost_per_case,0.75),'%4.1f') ')'];

Cost_per_Case{5} = [num2str(icdf(pd_4_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_4_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_4_cost_per_case,0.75),'%4.1f') ')'];

Cost_per_Case{6} = [num2str(icdf(pd_5_cost_per_case,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_5_cost_per_case,0.25),'%4.1f') char(8211) num2str(icdf(pd_5_cost_per_case,0.75),'%4.1f') ')'];


Output_Table=table(Year,Cases,Hospitalizations,Cost,Cost_per_Case);

close all;
figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

Current_Case_Count=1723;
MDL=icdf(pd_cases,[0.5 0.25 0.75]);
s1=subplot("Position",[0.075 0.15 0.87 0.8]);
cases=linspace(0,7500,10001);
plot(cases,pdf(pd_cases,cases),'color','k','LineWidth',2); hold on;
ym=max(pdf(pd_cases,cases));


p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none');
p1=plot([Current_Case_Count Current_Case_Count],[0 ym.*1.05],'r-.','LineWidth',2);
% p2=plot([Rough_Estimate Rough_Estimate],[0 ym.*1.05],'m-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);

% legend([p1 p2 p3 p4],{'Current count','Approximate true count','Model median','Model IQR'})
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cases(1):750:10^4],'Fontsize',16)
xlabel(['Annual measles cases'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([0,cases(end)])
ylim([0 ym.*1.05])
print(gcf,['HELP_Incidence_National.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.075 0.15 0.87 0.8]);
hospital=linspace(0,800,10001);
plot(hospital,pdf(pd_hospital,hospital),'color','k','LineWidth',2); hold on;
set(gca,'LineWidth',2,'Tickdir','out','XTick',[hospital(1):50:hospital(end)],'Fontsize',16)
xlabel(['Annual measles hospitalizations'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([hospital(1),hospital(end)])
print(gcf,['HELP_Hospital_National.png'],'-dpng','-r300');

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.095 0.15 0.87 0.8]);
cost=linspace(0,300,10001);
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
% 
