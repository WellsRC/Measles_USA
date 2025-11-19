clear;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Baseline Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
[pd_cases]=National_Outcome_Distribution(0,'Sample_2024',0);

close all;
figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

Current_Case_Count=285;
MDL=icdf(pd_cases,[0.5 0.25 0.75]);
s1=subplot("Position",[0.075 0.15 0.87 0.8]);
cases=linspace(0,800,10001);
plot(cases,pdf(pd_cases,cases),'color','k','LineWidth',2); hold on;
ym=max(pdf(pd_cases,cases));


% p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none');
p1=plot([Current_Case_Count Current_Case_Count],[0 ym.*1.05],'r-.','LineWidth',2);
% p2=plot([Rough_Estimate Rough_Estimate],[0 ym.*1.05],'m-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);

% legend([p1 p2 p3 p4],{'Current count','Approximate true count','Model median','Model IQR'})
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cases(1):100:800],'Fontsize',16)
xlabel(['Annual measles cases'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([0,cases(end)])
ylim([0 ym.*1.05])
print(gcf,['Validate_2024.png'],'-dpng','-r300');


[pd_cases]=National_Outcome_Distribution(0,'Sample_2023',0);
figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

Current_Case_Count=59;
MDL=icdf(pd_cases,[0.5 0.25 0.75]);
s1=subplot("Position",[0.105 0.15 0.86 0.8]);
cases=linspace(0,250,10001);
plot(cases,pdf(pd_cases,cases),'color','k','LineWidth',2); hold on;
ym=max(pdf(pd_cases,cases));


% p4=patch([MDL(2) MDL(2) MDL(3) MDL(3)],[0 ym.*1.05 ym.*1.05 0],'k','FaceAlpha',0.3,'LineStyle','none');
p1=plot([Current_Case_Count Current_Case_Count],[0 ym.*1.05],'r-.','LineWidth',2);
% p2=plot([Rough_Estimate Rough_Estimate],[0 ym.*1.05],'m-.','LineWidth',2);
p3=plot([MDL(1) MDL(1)],[0 ym.*1.05],'k-.','LineWidth',2);

% legend([p1 p2 p3 p4],{'Current count','Approximate true count','Model median','Model IQR'})
set(gca,'LineWidth',2,'Tickdir','out','XTick',[cases(1):25:250],'Fontsize',16)
xlabel(['Annual measles cases'],'FontSize',18)
ylabel('Probability density','FontSize',18)
box off;
xlim([0,cases(end)])
ylim([0 ym.*1.05])
print(gcf,['Validate_2023.png'],'-dpng','-r300');