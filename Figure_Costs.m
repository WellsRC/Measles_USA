function Figure_Costs(Scenario_Plot,Age_0_to_6)

close all;
        
    if(strcmp(Scenario_Plot,'Baseline'))
        Xlabel_Add_TXT={'Baseline'};
    elseif(strcmp(Scenario_Plot,'Sample'))
        Xlabel_Add_TXT={'Generalized'};
    else
        Xlabel_Add_TXT={num2str(Scenario_Plot(end-3:end))};
    end

    
    if(~Age_0_to_6)
        Age_Reduction=[true(1,1) false(1,4)];
    else
        Age_Reduction=[true(1,2) false(1,3)];
    end
    NS=2.5.*10^3;
    Colors={'#fecc5c';
'#fd8d3c';
'#e31a1c';};
    CC=hex2rgb(Colors);

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Baseline Calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    National_Reduction=0;
    [~,~,~,~,pd_pro_loss_baseline,pd_med_cost_baseline,pd_test_vac_cost_baseline,pd_ct_cost_baseline,pd_pro_loss_per_case_baseline,pd_med_cost_per_case_baseline,pd_test_vac_cost_per_case_baseline,pd_ct_cost_per_case_baseline]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);
    
    National_Reduction=0.01;
    [~,~,~,~,pd_pro_loss_reduction_1,pd_med_cost_reduction_1,pd_test_vac_cost_reduction_1,pd_ct_cost_reduction_1,pd_pro_loss_per_case_reduction_1,pd_med_cost_per_case_reduction_1,pd_test_vac_cost_per_case_reduction_1,pd_ct_cost_per_case_reduction_1]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);

    National_Reduction=0.025;
    [~,~,~,~,pd_pro_loss_reduction_2_5,pd_med_cost_reduction_2_5,pd_test_vac_cost_reduction_2_5,pd_ct_cost_reduction_2_5,pd_pro_loss_per_case_reduction_2_5,pd_med_cost_per_case_reduction_2_5,pd_test_vac_cost_per_case_reduction_2_5,pd_ct_cost_per_case_reduction_2_5]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);

    National_Reduction=0.05;
    [~,~,~,~,pd_pro_loss_reduction_5,pd_med_cost_reduction_5,pd_test_vac_cost_reduction_5,pd_ct_cost_reduction_5,pd_pro_loss_per_case_reduction_5,pd_med_cost_per_case_reduction_5,pd_test_vac_cost_per_case_reduction_5,pd_ct_cost_per_case_reduction_5]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);
 
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    s1=subplot("Position",[0.075 0.605 0.4 0.35]);

    [pro_loss,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_pro_loss_baseline,pd_pro_loss_reduction_1,pd_pro_loss_reduction_2_5,pd_pro_loss_reduction_5);   

    plot(pro_loss./10^6,pdf(pd_pro_loss_baseline,pro_loss),'k','LineWidth',2); hold on;
    plot(pro_loss./10^6,pdf(pd_pro_loss_reduction_1,pro_loss),'color',CC(1,:),'LineWidth',2);
    plot(pro_loss./10^6,pdf(pd_pro_loss_reduction_2_5,pro_loss),'color',CC(2,:),'LineWidth',2);
    plot(pro_loss./10^6,pdf(pd_pro_loss_reduction_5,pro_loss),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Productivity losses (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,max_x]./10^6)
    xtickformat('$%,.0f')

    myAxes=findobj(s1,'Type','Axes');
    exportgraphics(myAxes,['National_Total_Productivity_Loss_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'A','Fontsize',40,'Units','Normalized')

    s2=subplot("Position",[0.575 0.605 0.4 0.35]);
    [hospital,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_med_cost_baseline,pd_med_cost_reduction_1,pd_med_cost_reduction_2_5,pd_med_cost_reduction_5);
    
    plot(hospital./10^6,pdf(pd_med_cost_baseline,hospital),'k','LineWidth',2); hold on;
    plot(hospital./10^6,pdf(pd_med_cost_reduction_1,hospital),'color',CC(1,:),'LineWidth',2);
    plot(hospital./10^6,pdf(pd_med_cost_reduction_2_5,hospital),'color',CC(2,:),'LineWidth',2);
    plot(hospital./10^6,pdf(pd_med_cost_reduction_5,hospital),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Direct medical costs (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,max_x]./10^6)
    xtickformat('$%,.0f')
    
    myAxes=findobj(s2,'Type','Axes');
    exportgraphics(myAxes,['National_Total_Direct_Medical_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'B','Fontsize',40,'Units','Normalized')

    s3=subplot("Position",[0.075 0.105 0.4 0.35]);
    [Test_Vac_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_test_vac_cost_baseline,pd_test_vac_cost_reduction_1,pd_test_vac_cost_reduction_2_5,pd_test_vac_cost_reduction_5);
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_baseline,Test_Vac_cost),'k','LineWidth',2); hold on;
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_reduction_1,Test_Vac_cost),'color',CC(1,:),'LineWidth',2);
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_reduction_2_5,Test_Vac_cost),'color',CC(2,:),'LineWidth',2);
    plot(Test_Vac_cost./10^6,pdf(pd_test_vac_cost_reduction_5,Test_Vac_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Cost of testing and vaccinatiing contacts (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x]./10^6)

    myAxes=findobj(s3,'Type','Axes');
    exportgraphics(myAxes,['National_Cost_Test_Vaccination_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'C','Fontsize',40,'Units','Normalized')

    s4=subplot("Position",[0.575 0.105 0.4 0.35]);
    [ct_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_ct_cost_baseline,pd_ct_cost_reduction_1,pd_ct_cost_reduction_2_5,pd_ct_cost_reduction_5);
    plot(ct_cost./10^6,pdf(pd_ct_cost_baseline,ct_cost),'k','LineWidth',2); hold on;
    plot(ct_cost./10^6,pdf(pd_ct_cost_reduction_1,ct_cost),'color',CC(1,:),'LineWidth',2);
    plot(ct_cost./10^6,pdf(pd_ct_cost_reduction_2_5,ct_cost),'color',CC(2,:),'LineWidth',2);
    plot(ct_cost./10^6,pdf(pd_ct_cost_reduction_5,ct_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^6,'Fontsize',16)
    xlabel(['Cost of contact tracing (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x]./10^6)

    myAxes=findobj(s4,'Type','Axes');
    exportgraphics(myAxes,['National_CT_Cost_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'D','Fontsize',40,'Units','Normalized')

    print(gcf,['National_Total_Cost_Meaures_' Xlabel_Add_TXT{:} '.png'],'-dpng','-r300');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Per Case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    s1=subplot("Position",[0.075 0.605 0.4 0.35]);

    [pro_loss,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_pro_loss_per_case_baseline,pd_pro_loss_per_case_reduction_1,pd_pro_loss_per_case_reduction_2_5,pd_pro_loss_per_case_reduction_5);   

    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_baseline,pro_loss),'k','LineWidth',2); hold on;
    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_reduction_1,pro_loss),'color',CC(1,:),'LineWidth',2);
    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_reduction_2_5,pro_loss),'color',CC(2,:),'LineWidth',2);
    plot(pro_loss./10^3,pdf(pd_pro_loss_per_case_reduction_5,pro_loss),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^3,'Fontsize',16)
    xlabel(['Productivity losses per case (Thousands)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,max_x]./10^3)
    xtickformat('$%,.0f')

    myAxes=findobj(s1,'Type','Axes');
    exportgraphics(myAxes,['National_Per_Case_Productivity_Loss_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'A','Fontsize',40,'Units','Normalized')

    s2=subplot("Position",[0.575 0.605 0.4 0.35]);
    [hospital,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_med_cost_per_case_baseline,pd_med_cost_per_case_reduction_1,pd_med_cost_per_case_reduction_2_5,pd_med_cost_per_case_reduction_5);
    
    plot(hospital,pdf(pd_med_cost_per_case_baseline,hospital),'k','LineWidth',2); hold on;
    plot(hospital,pdf(pd_med_cost_per_case_reduction_1,hospital),'color',CC(1,:),'LineWidth',2);
    plot(hospital,pdf(pd_med_cost_per_case_reduction_2_5,hospital),'color',CC(2,:),'LineWidth',2);
    plot(hospital,pdf(pd_med_cost_per_case_reduction_5,hospital),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Direct medical costs per case'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,max_x])
    xtickformat('$%,.0f')
    
    myAxes=findobj(s2,'Type','Axes');
    exportgraphics(myAxes,['National_Per_Case_Direct_Medical_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'B','Fontsize',40,'Units','Normalized')

    s3=subplot("Position",[0.075 0.105 0.4 0.35]);
    [Test_Vac_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_test_vac_cost_per_case_baseline,pd_test_vac_cost_per_case_reduction_1,pd_test_vac_cost_per_case_reduction_2_5,pd_test_vac_cost_per_case_reduction_5);
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_baseline,Test_Vac_cost),'k','LineWidth',2); hold on;
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_reduction_1,Test_Vac_cost),'color',CC(1,:),'LineWidth',2);
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_reduction_2_5,Test_Vac_cost),'color',CC(2,:),'LineWidth',2);
    plot(Test_Vac_cost,pdf(pd_test_vac_cost_per_case_reduction_5,Test_Vac_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:100:max_x],'Fontsize',16)
    xlabel(['Cost of testing and vaccinatiing contacts per case'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x])

    myAxes=findobj(s3,'Type','Axes');
    exportgraphics(myAxes,['National_Per_Case_Cost_Test_Vaccination_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'C','Fontsize',40,'Units','Normalized')

    s4=subplot("Position",[0.575 0.105 0.4 0.35]);
    [ct_cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_ct_cost_per_case_baseline,pd_ct_cost_per_case_reduction_1,pd_ct_cost_per_case_reduction_2_5,pd_ct_cost_per_case_reduction_5);
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_baseline,ct_cost),'k','LineWidth',2); hold on;
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_reduction_1,ct_cost),'color',CC(1,:),'LineWidth',2);
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_reduction_2_5,ct_cost),'color',CC(2,:),'LineWidth',2);
    plot(ct_cost./10^3,pdf(pd_ct_cost_per_case_reduction_5,ct_cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x]./10^3,'Fontsize',16)
    xlabel(['Cost of contact tracing per case (Thousands)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.1f')
    xlim([min_x,max_x]./10^3)

    myAxes=findobj(s4,'Type','Axes');
    exportgraphics(myAxes,['National_CT_Cost_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'D','Fontsize',40,'Units','Normalized')

    print(gcf,['National_Per_Case_Cost_Meaures_' Xlabel_Add_TXT{:} '.png'],'-dpng','-r300');
end
