clear;
close all;
clc;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Baseline Calculations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 

Measure={'Cases';'Hospitalizations';'Cost';'Outbreak respons';'Direct Medical';'Direct Medical: Uninsured';'Direct Medical: Publicly insured';'Direct Medical: Privately insured';'Productivity loses';'Cost per case'};
[pd_cases_2025,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_med_cost_uninsured,pd_med_cost_public,pd_med_cost_private,~,~,pd_outbreak_response_cost,pd_severe_disease]=National_Outcome_Distribution(0,'Baseline',0);

[num2str(icdf(pd_severe_disease,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_severe_disease,0.25),'%5.0f') char(8211) num2str(icdf(pd_severe_disease,0.75),'%5.0f') ')']

pd_cases=pd_cases_2025;
Value_2025=cell(10,1);
  Value_2025{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];

    Value_2025{2} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];

    Value_2025{3} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ') million'];

    Value_2025{4} = [num2str(icdf(pd_outbreak_response_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_outbreak_response_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_outbreak_response_cost,0.75)./10^6,'%4.1f') ') million'];

    Value_2025{5} = [num2str(icdf(pd_med_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_med_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_med_cost,0.75)./10^6,'%4.1f') ') million'];

    Value_2025{6} = [num2str(icdf(pd_med_cost_uninsured,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_uninsured,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_uninsured,0.75)./10^6,'%4.2f') ') million'];
    Value_2025{7} = [num2str(icdf(pd_med_cost_public,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_public,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_public,0.75)./10^6,'%4.2f') ') million'];
    Value_2025{8} = [num2str(icdf(pd_med_cost_private,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_private,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_private,0.75)./10^6,'%4.2f') ') million'];

    Value_2025{9} = [num2str(icdf(pd_pro_loss,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_pro_loss,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_pro_loss,0.75)./10^6,'%4.1f') ') million'];

    Value_2025{10} = [num2str(10^3.*icdf(pd_cost_per_case,0.5),'%5.0f') ' (IQR:' num2str(10^3.*icdf(pd_cost_per_case,0.25),'%5.0f') char(8211) num2str(10^3.*icdf(pd_cost_per_case,0.75),'%5.0f') ')'];

[pd_cases_2024,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_med_cost_uninsured,pd_med_cost_public,pd_med_cost_private,~,~,pd_outbreak_response_cost]=National_Outcome_Distribution(0,'Sample_2024',0);

pd_cases=pd_cases_2024;
Value_2024=cell(10,1);
  Value_2024{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];

    Value_2024{2} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];

    Value_2024{3} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ') million'];

    Value_2024{4} = [num2str(icdf(pd_outbreak_response_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_outbreak_response_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_outbreak_response_cost,0.75)./10^6,'%4.1f') ') million'];

    Value_2024{5} = [num2str(icdf(pd_med_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_med_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_med_cost,0.75)./10^6,'%4.1f') ') million'];

     Value_2024{6} = [num2str(icdf(pd_med_cost_uninsured,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_uninsured,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_uninsured,0.75)./10^6,'%4.2f') ') million'];
    Value_2024{7} = [num2str(icdf(pd_med_cost_public,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_public,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_public,0.75)./10^6,'%4.2f') ') million'];
    Value_2024{8} = [num2str(icdf(pd_med_cost_private,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_private,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_private,0.75)./10^6,'%4.2f') ') million'];

    Value_2024{9} = [num2str(icdf(pd_pro_loss,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_pro_loss,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_pro_loss,0.75)./10^6,'%4.1f') ') million'];

    Value_2024{10} = [num2str(10^3.*icdf(pd_cost_per_case,0.5),'%5.0f') ' (IQR:' num2str(10^3.*icdf(pd_cost_per_case,0.25),'%5.0f') char(8211) num2str(10^3.*icdf(pd_cost_per_case,0.75),'%5.0f') ')'];


[pd_cases_2023,pd_hospital,pd_cost,pd_cost_per_case,pd_pro_loss,pd_med_cost,pd_med_cost_uninsured,pd_med_cost_public,pd_med_cost_private,~,~,pd_outbreak_response_cost]=National_Outcome_Distribution(0,'Sample_2023',0);

pd_cases=pd_cases_2023;
Value_2023=cell(10,1);
  Value_2023{1} = [num2str(icdf(pd_cases,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_cases,0.25),'%5.0f') char(8211) num2str(icdf(pd_cases,0.75),'%5.0f') ')'];
    
    Value_2023{2} = [num2str(icdf(pd_hospital,0.5),'%5.0f') ' (IQR:' num2str(icdf(pd_hospital,0.25),'%5.0f') char(8211) num2str(icdf(pd_hospital,0.75),'%5.0f') ')'];
    
    Value_2023{3} = [num2str(icdf(pd_cost,0.5),'%4.1f') ' (IQR:' num2str(icdf(pd_cost,0.25),'%4.1f') char(8211) num2str(icdf(pd_cost,0.75),'%4.1f') ') million'];
    
    Value_2023{4} = [num2str(icdf(pd_outbreak_response_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_outbreak_response_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_outbreak_response_cost,0.75)./10^6,'%4.1f') ') million'];
    
    Value_2023{5} = [num2str(icdf(pd_med_cost,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_med_cost,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_med_cost,0.75)./10^6,'%4.1f') ') million'];

     Value_2023{6} = [num2str(icdf(pd_med_cost_uninsured,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_uninsured,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_uninsured,0.75)./10^6,'%4.2f') ') million'];
    Value_2023{7} = [num2str(icdf(pd_med_cost_public,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_public,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_public,0.75)./10^6,'%4.2f') ') million'];
    Value_2023{8} = [num2str(icdf(pd_med_cost_private,0.5)./10^6,'%4.2f') ' (IQR:' num2str(icdf(pd_med_cost_private,0.25)./10^6,'%4.2f') char(8211) num2str(icdf(pd_med_cost_private,0.75)./10^6,'%4.2f') ') million'];


    Value_2023{9} = [num2str(icdf(pd_pro_loss,0.5)./10^6,'%4.1f') ' (IQR:' num2str(icdf(pd_pro_loss,0.25)./10^6,'%4.1f') char(8211) num2str(icdf(pd_pro_loss,0.75)./10^6,'%4.1f') ') million'];

    Value_2023{10} = [num2str(10^3.*icdf(pd_cost_per_case,0.5),'%5.0f') ' (IQR:' num2str(10^3.*icdf(pd_cost_per_case,0.25),'%5.0f') char(8211) num2str(10^3.*icdf(pd_cost_per_case,0.75),'%5.0f') ')'];

 Output_Table=table(Measure,Value_2025,Value_2024,Value_2023);


close all;
f=figure('units','normalized','outerposition',[0.05 0.05 0.8 1]);
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

theme(f,'light');
print(gcf,['Figure_S2.png'],'-dpng','-r300');
