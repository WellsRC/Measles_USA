clear;
close all;
clc;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Baseline Calculations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
[pd_cases_2025]=National_Outcome_Distribution(0,'Baseline',0);
[pd_cases_2024]=National_Outcome_Distribution(0,'Sample_2024',0);
[pd_cases_2023]=National_Outcome_Distribution(0,'Sample_2023',0);


close all;
figure('units','normalized','outerposition',[0.05 0.05 0.8 1]);
s1=subplot("Position",[0.075 0.58 0.87 0.39]);

Case_Count_2025=1753;
MDL=icdf(pd_cases_2025,[0.5 0.25 0.75]);
cases=linspace(0,3750,10001);
ym=max(pdf(pd_cases_2025,cases));
p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none'); hold on;
plot(cases,pdf(pd_cases_2025,cases),'color','k','LineWidth',2);
p1=plot([Case_Count_2025 Case_Count_2025],[0 ym.*1.05],'r-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);
ylim([0 ym.*1.05])
xlim([0 cases(end)])
set(gca,'LineWidth',2,'Tickdir','out','XTick',[0:250:3750],'Fontsize',16)
xlabel(['Annual measles cases (2025)'],'FontSize',18)
ylabel('Probability density','FontSize',18)

text(-0.06,1,'A','FontSize',34,'Units','normalized');

s2=subplot("Position",[0.075 0.08 0.4 0.4]);

Case_Count_2024=285;
MDL=icdf(pd_cases_2024,[0.5 0.25 0.75]);
cases=linspace(0,650,10001);
ym=max(pdf(pd_cases_2024,cases));
p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none'); hold on;
plot(cases,pdf(pd_cases_2024,cases),'color','k','LineWidth',2);
p1=plot([Case_Count_2024 Case_Count_2024],[0 ym.*1.05],'r-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);
ylim([0 ym.*1.05])
xlim([0 cases(end)])

set(gca,'LineWidth',2,'Tickdir','out','XTick',[0:75:650],'Fontsize',16)
xlabel(['Annual measles cases (2024)'],'FontSize',18)
ylabel('Probability density','FontSize',18)

text(-0.13,1,'B','FontSize',34,'Units','normalized');

s3=subplot("Position",[0.55 0.08 0.4 0.4]);

Case_Count_2023=59;
MDL=icdf(pd_cases_2023,[0.5 0.25 0.75]);
cases=linspace(0,200,10001);
ym=max(pdf(pd_cases_2023,cases));
p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none'); hold on;
plot(cases,pdf(pd_cases_2023,cases),'color','k','LineWidth',2);
p1=plot([Case_Count_2023 Case_Count_2023],[0 ym.*1.05],'r-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);
ylim([0 ym.*1.05])
xlim([0 cases(end)])
set(gca,'LineWidth',2,'Tickdir','out','XTick',[0:20:200],'Fontsize',16)
xlabel(['Annual measles cases (2023)'],'FontSize',18)
ylabel('Probability density','FontSize',18)
text(-0.1711,1,'C','FontSize',34,'Units','normalized');
print(gcf,['Figure_2.png'],'-dpng','-r300');
