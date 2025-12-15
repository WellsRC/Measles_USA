clear;
clc;
close all;
f=figure('units','normalized','outerposition',[0.05 0.05 0.6 1]);
CV=zeros(length(25:25:200),1);
for imprt=25:25:200
    subplot('Position',[0.125,0.58,0.85,0.395]);
    load(['Monte_Carlo_Run_Scenario_' num2str(imprt) '_National_Reduction=0_Year=0.mat']);

    Samp_B=sum(Total_Cases_County,1);
    
    CV(imprt./25)=std(Samp_B)./mean(Samp_B);
    
    Cases_Baseline=prctile(sum(Total_Cases_County,1),[50 25 75]);
    patch(imprt+[-11 -11 11 11],[Cases_Baseline(2) Cases_Baseline(3) Cases_Baseline(3) Cases_Baseline(2)],'k','LineStyle','none','FaceAlpha',0.3); hold on
    plot(imprt+[-11 11],[Cases_Baseline(1) Cases_Baseline(1)],'k','LineWidth',2);


    load(['Monte_Carlo_Run_Scenario_' num2str(imprt) '_National_Reduction=0.25_Year=1.mat']);
    Samp_R1=sum(Total_Cases_County,1);
    
    load(['Monte_Carlo_Run_Scenario_' num2str(imprt) '_National_Reduction=0.5_Year=1.mat']);
    Samp_R2=sum(Total_Cases_County,1);
    
    load(['Monte_Carlo_Run_Scenario_' num2str(imprt) '_National_Reduction=0.75_Year=1.mat']);
    Samp_R3=sum(Total_Cases_County,1);
    
    load(['Monte_Carlo_Run_Scenario_' num2str(imprt) '_National_Reduction=1_Year=1.mat']);
    Samp_R4=sum(Total_Cases_County,1);
    
    M1=median(Samp_R1-Samp_B);
    M2=median(Samp_R2-Samp_B);
    M3=median(Samp_R3-Samp_B);
    M4=median(Samp_R4-Samp_B);
    subplot('Position',[0.125,0.078,0.85,0.395]);
    xx=linspace(imprt-11,imprt+11,5);
    dx=0.4.*(xx(2)-xx(1));
    xx=(xx(2:end)+xx(1:end-1))./2;
    p1=patch(xx(1)+[-dx -dx dx dx],[0 M1 M1 0],hex2rgb('#ffffff'),'LineWidth',2,'EdgeColor',hex2rgb('#67000d')); hold on
    p2=patch(xx(2)+[-dx -dx dx dx],[0 M2 M2 0],hex2rgb('#fc9272'),'LineWidth',2,'EdgeColor',hex2rgb('#fc9272')); hold on
    p3=patch(xx(3)+[-dx -dx dx dx],[0 M3 M3 0],hex2rgb('#cb181d'),'LineWidth',2,'EdgeColor',hex2rgb('#cb181d')); hold on
    p4=patch(xx(4)+[-dx -dx dx dx],[0 M4 M4 0],hex2rgb('#67000d'),'LineWidth',2,'EdgeColor',hex2rgb('#67000d')); hold on

end

subplot('Position',[0.125,0.58,0.85,0.395]);
set(gca,'XTick',[25:25:200],'tickdir','out','lineWidth',2,'yscale','log','FontSize',16)
xlabel('Annual importation events','FontSize',18)
ylabel('Annual measles incidence','FontSize',18)
ylim([10 10^5])
text(-0.12,1,'A','Fontsize',32,'Units','normalized')
xlim([12.5 212.5])
subplot('Position',[0.125,0.078,0.85,0.395]);
set(gca,'XTick',[25:25:200],'tickdir','out','lineWidth',2,'yscale','linear','FontSize',16)
xlabel('Annual importation events','FontSize',18)
ylabel('Additional annual measles incidence','FontSize',18)
legend([p1 p2 p3 p4],{'0.25%','0.50%','0.75%','1.00%'},'FontSize',16,'Location','northwest')
text(-0.12,1,'B','Fontsize',32,'Units','normalized')
xlim([12.5 212.5])

theme(f, "light"); 
print(gcf,['Figure_2.png'],'-dpng','-r300');
print(gcf,['Figure_2.tif'],'-dtiff','-r300');