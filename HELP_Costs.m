clear;
clc;
close all;

if(~Age_0_to_6)
    Age_Reduction=[true(1,1) false(1,4)];
else
    Age_Reduction=[true(1,2) false(1,3)];
end

Colors={'#1E434C';
'#C99E10';
'#8D230F';};
CC=hex2rgb(Colors);


    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Baseline Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HELP_Reduction=0.01;
    [~,~,~,~,pd_pro_loss_2023,pd_med_cost_2023,pd_test_vac_cost_2023,pd_ct_cost_2023,pd_pro_loss_per_case_2023,pd_med_cost_per_case_2023,pd_test_vac_cost_per_case_2023,pd_ct_cost_per_case_2023]=HELP_Outcome_Distribution(HELP_Reduction,'Sample_2023',Age_0_to_6);
    [~,~,~,~,pd_pro_loss_2024,pd_med_cost_2024,pd_test_vac_cost_2024,pd_ct_cost_2024,pd_pro_loss_per_case_2024,pd_med_cost_per_case_2024,pd_test_vac_cost_per_case_2024,pd_ct_cost_per_case_2024]=HELP_Outcome_Distribution(HELP_Reduction,'Sample_2024',Age_0_to_6);
    [~,~,~,~,pd_pro_loss_2025,pd_med_cost_2025,pd_test_vac_cost_2025,pd_ct_cost_2025,pd_pro_loss_per_case_2025,pd_med_cost_per_case_2025,pd_test_vac_cost_per_case_2025,pd_ct_cost_per_case_2025]=HELP_Outcome_Distribution(HELP_Reduction,'Sample_2025',Age_0_to_6);
 
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    s1=subplot("Position",[0.075 0.605 0.4 0.35]);

    [pro_loss,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_pro_loss_baseline,pd_pro_loss_2023,pd_pro_loss_2024,pd_pro_loss_2025);   

    plot(pro_loss./10^6,pdf(pd_pro_loss_2023,pro_loss),'color',CC(1,:),'LineWidth',2); hold on;
    plot(pro_loss./10^6,pdf(pd_pro_loss_2024,pro_loss),'color',CC(2,:),'LineWidth',2);
    plot(pro_loss./10^6,pdf(pd_pro_loss_2025,pro_loss),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Productivity losses (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xlim([min_x,max_x]./10^6)
    xtickformat('$%,.0f')

    myAxes=findobj(s1,'Type','Axes');
    exportgraphics(myAxes,['HELP_Total_Productivity_Loss.png'],'Resolution',300);

    text(-0.15,1.05,'A','Fontsize',40,'Units','Normalized')

    s2=subplot("Position",[0.575 0.605 0.4 0.35]);
    [hospital,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_med_cost_baseline,pd_med_cost_2023,pd_med_cost_2024,pd_med_cost_2025);
    
    plot(hospital./10^6,pdf(pd_med_cost_2023,hospital),'color',CC(1,:),'LineWidth',2); hold on;
    plot(hospital./10^6,pdf(pd_med_cost_2024,hospital),'color',CC(2,:),'LineWidth',2);
    plot(hospital./10^6,pdf(pd_med_cost_2025,hospital),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Direct medical costs (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xlim([min_x,max_x]./10^6)
    xtickformat('$%,.0f')
    
    myAxes=findobj(s2,'Type','Axes');
    exportgraphics(myAxes,['HELP_Total_Direct_Medical.png'],'Resolution',300);

    text(-0.15,1.05,'B','Fontsize',40,'Units','Normalized')

    s3=subplot("Position",[0.075 0.105 0.4 0.35]);
    [Test_Vac_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_test_vac_cost_baseline,pd_test_vac_cost_2023,pd_test_vac_cost_2024,pd_test_vac_cost_2025);
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_2023,Test_Vac_cost),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_2024,Test_Vac_cost),'color',CC(2,:),'LineWidth',2);
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_2025,Test_Vac_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Cost of testing and vaccinatiing contacts (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x]./10^6)

    myAxes=findobj(s3,'Type','Axes');
    exportgraphics(myAxes,['HELP_Cost_Test_Vaccination.png'],'Resolution',300);

    text(-0.15,1.05,'C','Fontsize',40,'Units','Normalized')

    s4=subplot("Position",[0.575 0.105 0.4 0.35]);
    [ct_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_ct_cost_baseline,pd_ct_cost_2023,pd_ct_cost_2024,pd_ct_cost_2025);
        plot(ct_cost./10^6,pdf(pd_ct_cost_2023,ct_cost),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(ct_cost./10^6,pdf(pd_ct_cost_2024,ct_cost),'color',CC(2,:),'LineWidth',2);
    plot(ct_cost./10^6,pdf(pd_ct_cost_2025,ct_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Cost of contact tracing (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x]./10^6)

    myAxes=findobj(s4,'Type','Axes');
    exportgraphics(myAxes,['HELP_CT_Cost.png'],'Resolution',300);

    text(-0.15,1.05,'D','Fontsize',40,'Units','Normalized')

    print(gcf,['HELP_Total_Cost_Meaures.png'],'-dpng','-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Per Case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    s1=subplot("Position",[0.075 0.605 0.4 0.35]);

    [pro_loss,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_pro_loss_per_case_baseline,pd_pro_loss_per_case_2023,pd_pro_loss_per_case_2024,pd_pro_loss_per_case_2025);   

    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_2023,pro_loss),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_2024,pro_loss),'color',CC(2,:),'LineWidth',2);
    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_2025,pro_loss),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^3,'Fontsize',16)
    xlabel(['Productivity losses per case (Thousands)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xlim([min_x,max_x]./10^3)
    xtickformat('$%,.0f')

    myAxes=findobj(s1,'Type','Axes');
    exportgraphics(myAxes,['HELP_Per_Case_Productivity_Loss.png'],'Resolution',300);

    text(-0.15,1.05,'A','Fontsize',40,'Units','Normalized')

    s2=subplot("Position",[0.575 0.605 0.4 0.35]);
    [hospital,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_med_cost_per_case_baseline,pd_med_cost_per_case_2023,pd_med_cost_per_case_2024,pd_med_cost_per_case_2025);
    
    plot(hospital,pdf(pd_med_cost_per_case_2023,hospital),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(hospital,pdf(pd_med_cost_per_case_2024,hospital),'color',CC(2,:),'LineWidth',2);
    plot(hospital,pdf(pd_med_cost_per_case_2025,hospital),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Direct medical costs per case'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xlim([min_x,max_x])
    xtickformat('$%,.0f')
    
    myAxes=findobj(s2,'Type','Axes');
    exportgraphics(myAxes,['HELP_Per_Case_Direct_Medical.png'],'Resolution',300);

    text(-0.15,1.05,'B','Fontsize',40,'Units','Normalized')

    s3=subplot("Position",[0.075 0.105 0.4 0.35]);
    [Test_Vac_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_test_vac_cost_per_case_baseline,pd_test_vac_cost_per_case_2023,pd_test_vac_cost_per_case_2024,pd_test_vac_cost_per_case_2025);
    
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_2023,Test_Vac_cost),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_2024,Test_Vac_cost),'color',CC(2,:),'LineWidth',2);
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_2025,Test_Vac_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:100:max_x],'Fontsize',16)
    xlabel(['Cost of testing and vaccinatiing contacts per case'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x])

    myAxes=findobj(s3,'Type','Axes');
    exportgraphics(myAxes,['HELP_Per_Case_Cost_Test_Vaccination.png'],'Resolution',300);

    text(-0.15,1.05,'C','Fontsize',40,'Units','Normalized')

    s4=subplot("Position",[0.575 0.105 0.4 0.35]);
    [ct_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_ct_cost_per_case_baseline,pd_ct_cost_per_case_2023,pd_ct_cost_per_case_2024,pd_ct_cost_per_case_2025);
    
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_2023,ct_cost),'color',CC(1,:),'LineWidth',2);  hold on;
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_2024,ct_cost),'color',CC(2,:),'LineWidth',2);
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_2025,ct_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^3,'Fontsize',16)
    xlabel(['Cost of contact tracing per case (Thousands)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    [hleg] = legend({'2023','2024','2025'},'FontSize',16,'Box','off');
    title(hleg,'Year')
    box off;
    xtickformat('$%,.1f')
    xlim([min_x,max_x]./10^3)

    myAxes=findobj(s4,'Type','Axes');
    exportgraphics(myAxes,['HELP_CT_Cost.png'],'Resolution',300);

    text(-0.15,1.05,'D','Fontsize',40,'Units','Normalized')

    print(gcf,['HELP_Per_Case_Cost_Meaures.png'],'-dpng','-r300');

