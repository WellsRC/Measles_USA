clear;
close all;

load('National_Reduction=0_Ages_0_to_6.mat')
X=1-County_Data_Vaccine_Reduction.Total_Immunity;
[cases,hospital,cost,cost_per_case]=County_Outcome_Central_Measure(0,'Sample_2025',true);
[r,p]=corr(X(:),cost_per_case(:),'Type','Spearman')

figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);

s1=subplot("Position",[0.095 0.15 0.87 0.8]);

scatter(100.*X(:),cost_per_case(:)/10^3,20,'k','filled')
xlabel('Susceptibility in the population','Fontsize',20)
ylabel('Cost per case (Thousands)','Fontsize',20)
ytickformat('$%,.0f')
xtickformat('percentage');
set(gca,'LineWidth',2,'tickdir','out','Fontsize',18)
box off
print(gcf,['Test_Immunity_vs_Cost_per_Case.png'],'-dpng','-r300');
