% clear;
% close all;
% clc;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Baseline Calculations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Age_0_to_6=true;
% 
% Reductions=[1 2 3 4 5 10 15 20 25];
% 
% figure('units','normalized','outerposition',[0.05 0.05 0.6 0.6]);
% 
% s1=subplot("Position",[0.095 0.15 0.87 0.8]);
% 
% for ii=1:length(Reductions)
%     [~,~,pd_cost]=National_Outcome_Comparison_Distribution(Reductions(ii)./100,'Sample_2025',Age_0_to_6);
% 
%     md=icdf(pd_cost,0.5);
%     lb2=icdf(pd_cost,0.025);
%     ub2=icdf(pd_cost,0.975);
%     lb1=icdf(pd_cost,0.25);
%     ub1=icdf(pd_cost,0.75);
% 
%     % patch(Reductions(ii)+[-0.25 -0.25 0.25 0.25],[lb2 ub2 ub2 lb2],'k','FaceAlpha',0.3,'LineStyle','none'); hold on
%     patch(Reductions(ii)+[-0.25 -0.25 0.25 0.25],[lb1 ub1 ub1 lb1],'k','FaceAlpha',0.3,'LineStyle','none'); hold on
%     plot(Reductions(ii)+[-0.25 0.25],[md md],'k','LineWidth',4);
% end

set(gca,'LineWidth',2,'Tickdir','out','XTick',[Reductions],'YTick',[0:25:175],'Fontsize',16)
xlabel(['Absoluat reduction in vaccine uptake (0' char(8211) '6 years of age)'],'FontSize',18)
xtickformat("percentage")
ylabel('Additional cost','FontSize',18)
ytickformat('$%,.0f')
box off;
xlim([0.5 25.5])
ylim([0 175])