function Output_Table=Figure_Incidence(Scenario_Plot,Age_0_to_6)
        
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
    [pd_baseline_cases,pd_baseline_hospital,pd_baseline_cost,pd_baseline_cost_per_case]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);
    
    National_Reduction=0.01;
    [pd_reduction_1_cases,pd_reduction_1_hospital,pd_reduction_1_cost,pd_reduction_1_cost_per_case]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);

    National_Reduction=0.025;
    [pd_reduction_2_5_cases,pd_reduction_2_5_hospital,pd_reduction_2_5_cost,pd_reduction_2_5_cost_per_case]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);

    National_Reduction=0.05;
    [pd_reduction_5_cases,pd_reduction_5_hospital,pd_reduction_5_cost,pd_reduction_5_cost_per_case]=National_Outcome_Distribution(National_Reduction,Scenario_Plot,Age_0_to_6);


    Scenario={'Baseline';'1% Reduction'; '2.5% Reduction'; '5% Reduction'};
    Cases=cell(length(Scenario),1);
    Hospitalizations=cell(length(Scenario),1);
    Cost=cell(length(Scenario),1);
    Cost_per_case=cell(length(Scenario),1);

    % Cases
    x0=icdf(pd_baseline_cases,0.5);
    lb=icdf(pd_baseline_cases,0.05);
    ub=icdf(pd_baseline_cases,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_baseline_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cases{1} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_baseline_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_baseline_cases,0.975),'%5.0f') ')'];

    x0=icdf(pd_reduction_1_cases,0.5);
    lb=icdf(pd_reduction_1_cases,0.05);
    ub=icdf(pd_reduction_1_cases,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_1_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cases{2} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_1_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_1_cases,0.975),'%5.0f') ')'];

    x0=icdf(pd_reduction_2_5_cases,0.5);
    lb=icdf(pd_reduction_2_5_cases,0.05);
    ub=icdf(pd_reduction_2_5_cases,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_2_5_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cases{3} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_2_5_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_2_5_cases,0.975),'%5.0f') ')'];

    x0=icdf(pd_reduction_5_cases,0.5);
    lb=icdf(pd_reduction_5_cases,0.05);
    ub=icdf(pd_reduction_5_cases,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_5_cases,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cases{4} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_5_cases,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_5_cases,0.975),'%5.0f') ')'];
    
    %Hospitalization
    x0=icdf(pd_baseline_hospital,0.5);
    lb=icdf(pd_baseline_hospital,0.05);
    ub=icdf(pd_baseline_hospital,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_baseline_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Hospitalizations{1} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_baseline_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_baseline_hospital,0.975),'%5.0f') ')'];
    
    x0=icdf(pd_reduction_1_hospital,0.5);
    lb=icdf(pd_reduction_1_hospital,0.05);
    ub=icdf(pd_reduction_1_hospital,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_1_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Hospitalizations{2} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_1_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_1_hospital,0.975),'%5.0f') ')'];
    
    x0=icdf(pd_reduction_2_5_hospital,0.5);
    lb=icdf(pd_reduction_2_5_hospital,0.05);
    ub=icdf(pd_reduction_2_5_hospital,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_2_5_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Hospitalizations{3} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_2_5_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_2_5_hospital,0.975),'%5.0f') ')'];
    
    x0=icdf(pd_reduction_5_hospital,0.5);
    lb=icdf(pd_reduction_5_hospital,0.05);
    ub=icdf(pd_reduction_5_hospital,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_5_hospital,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Hospitalizations{4} = [num2str(temp_X,'%5.0f') '(95% CrI:' num2str(icdf(pd_reduction_5_hospital,0.025),'%5.0f') char(8211) num2str(icdf(pd_reduction_5_hospital,0.975),'%5.0f') ')'];

    % Costs
    x0=icdf(pd_baseline_cost,0.5);
    lb=icdf(pd_baseline_cost,0.05);
    ub=icdf(pd_baseline_cost,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_baseline_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost{1} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_baseline_cost,0.025),'%4.2f') char(8211) num2str(icdf(pd_baseline_cost,0.975),'%4.2f') ')'];

    x0=icdf(pd_reduction_1_cost,0.5);
    lb=icdf(pd_reduction_1_cost,0.05);
    ub=icdf(pd_reduction_1_cost,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_1_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost{2} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_reduction_1_cost,0.025),'%4.2f') char(8211) num2str(icdf(pd_reduction_1_cost,0.975),'%4.2f') ')'];

    x0=icdf(pd_reduction_2_5_cost,0.5);
    lb=icdf(pd_reduction_2_5_cost,0.05);
    ub=icdf(pd_reduction_2_5_cost,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_2_5_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost{3} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_reduction_2_5_cost,0.025),'%4.2f') char(8211) num2str(icdf(pd_reduction_2_5_cost,0.975),'%4.2f') ')'];

    x0=icdf(pd_reduction_5_cost,0.5);
    lb=icdf(pd_reduction_5_cost,0.05);
    ub=icdf(pd_reduction_5_cost,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_5_cost,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost{4} = [num2str(temp_X,'%4.2f') '(95% CrI:' num2str(icdf(pd_reduction_5_cost,0.025),'%4.2f') char(8211) num2str(icdf(pd_reduction_5_cost,0.975),'%4.2f') ')'];


    % Costs per case
    x0=icdf(pd_baseline_cost_per_case,0.5);
    lb=icdf(pd_baseline_cost_per_case,0.05);
    ub=icdf(pd_baseline_cost_per_case,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_baseline_cost_per_case,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost_per_case{1} = [num2str(10.^3.*temp_X,'%5.0f') '(95% CrI:' num2str(10^3.*icdf(pd_baseline_cost_per_case,0.025),'%5.0f') char(8211) num2str(10^3.*icdf(pd_baseline_cost_per_case,0.975),'%5.0f') ')'];

    x0=icdf(pd_reduction_1_cost_per_case,0.5);
    lb=icdf(pd_reduction_1_cost_per_case,0.05);
    ub=icdf(pd_reduction_1_cost_per_case,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_1_cost_per_case,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost_per_case{2} = [num2str(10.^3.*temp_X,'%5.0f') '(95% CrI:' num2str(10^3.*icdf(pd_reduction_1_cost_per_case,0.025),'%5.0f') char(8211) num2str(10^3.*icdf(pd_reduction_1_cost_per_case,0.975),'%5.0f') ')'];
    
    x0=icdf(pd_reduction_2_5_cost_per_case,0.5);
    lb=icdf(pd_reduction_2_5_cost_per_case,0.05);
    ub=icdf(pd_reduction_2_5_cost_per_case,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_2_5_cost_per_case,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost_per_case{3} = [num2str(10.^3.*temp_X,'%5.0f') '(95% CrI:' num2str(10^3.*icdf(pd_reduction_2_5_cost_per_case,0.025),'%5.0f') char(8211) num2str(10^3.*icdf(pd_reduction_2_5_cost_per_case,0.975),'%5.0f') ')'];

    x0=icdf(pd_reduction_5_cost_per_case,0.5);
    lb=icdf(pd_reduction_5_cost_per_case,0.05);
    ub=icdf(pd_reduction_5_cost_per_case,0.95);
    temp_X=10.^fmincon(@(z)-log(pdf(pd_reduction_5_cost_per_case,10.^z)),log10(x0),[],[],[],[],log10(lb),log10(ub));
    Cost_per_case{4} = [num2str(10.^3.*temp_X,'%5.0f') '(95% CrI:' num2str(10^3.*icdf(pd_reduction_5_cost_per_case,0.025),'%5.0f') char(8211) num2str(10^3.*icdf(pd_reduction_5_cost_per_case,0.975),'%5.0f') ')'];

    Output_Table=table(Scenario,Cases,Hospitalizations,Cost,Cost_per_case);

    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    s1=subplot("Position",[0.075 0.605 0.4 0.35]);

    [cases,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_baseline_cases,pd_reduction_1_cases,pd_reduction_2_5_cases,pd_reduction_5_cases);   

    plot(cases,pdf(pd_baseline_cases,cases),'k','LineWidth',2); hold on;
    plot(cases,pdf(pd_reduction_1_cases,cases),'color',CC(1,:),'LineWidth',2);
    plot(cases,pdf(pd_reduction_2_5_cases,cases),'color',CC(2,:),'LineWidth',2);
    plot(cases,pdf(pd_reduction_5_cases,cases),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Annual measles cases'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,2500])
    
    myAxes=findobj(s1,'Type','Axes');
    exportgraphics(myAxes,['National_Incidence_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'A','Fontsize',40,'Units','Normalized')

    s2=subplot("Position",[0.575 0.605 0.4 0.35]);
    [hospital,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_baseline_hospital,pd_reduction_1_hospital,pd_reduction_2_5_hospital,pd_reduction_5_hospital);
    
    plot(hospital,pdf(pd_baseline_hospital,hospital),'k','LineWidth',2); hold on;
    plot(hospital,pdf(pd_reduction_1_hospital,hospital),'color',CC(1,:),'LineWidth',2);
    plot(hospital,pdf(pd_reduction_2_5_hospital,hospital),'color',CC(2,:),'LineWidth',2);
    plot(hospital,pdf(pd_reduction_5_hospital,hospital),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Annual measles hospitalizations'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xlim([min_x,700])
    
    myAxes=findobj(s2,'Type','Axes');
    exportgraphics(myAxes,['National_Hospital_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'B','Fontsize',40,'Units','Normalized')

    s3=subplot("Position",[0.075 0.105 0.4 0.35]);
    [cost,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_baseline_cost,pd_reduction_1_cost,pd_reduction_2_5_cost,pd_reduction_5_cost);
    plot(cost,pdf(pd_baseline_cost,cost),'k','LineWidth',2); hold on;
    plot(cost,pdf(pd_reduction_1_cost,cost),'color',CC(1,:),'LineWidth',2);
    plot(cost,pdf(pd_reduction_2_5_cost,cost),'color',CC(2,:),'LineWidth',2);
    plot(cost,pdf(pd_reduction_5_cost,cost),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Total costs (Millions)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,max_x])

    myAxes=findobj(s3,'Type','Axes');
    exportgraphics(myAxes,['National_Cost_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'C','Fontsize',40,'Units','Normalized')

    s4=subplot("Position",[0.575 0.105 0.4 0.35]);
    [cost_per_case,dx,min_x,max_x]=Bounds_for_Figure_Incidence(pd_baseline_cost_per_case,pd_reduction_1_cost_per_case,pd_reduction_2_5_cost_per_case,pd_reduction_5_cost_per_case);
    plot(cost_per_case,pdf(pd_baseline_cost_per_case,cost_per_case),'k','LineWidth',2); hold on;
    plot(cost_per_case,pdf(pd_reduction_1_cost_per_case,cost_per_case),'color',CC(1,:),'LineWidth',2);
    plot(cost_per_case,pdf(pd_reduction_2_5_cost_per_case,cost_per_case),'color',CC(2,:),'LineWidth',2);
    plot(cost_per_case,pdf(pd_reduction_5_cost_per_case,cost_per_case),'color',CC(3,:),'LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[min_x:dx:max_x],'Fontsize',16)
    xlabel(['Cost per case  (Thousands)'],'FontSize',18)
    ylabel('Probability density','FontSize',18)
    legend({'Baseline','1% reduction','2.5% reduction','5% reduction'},'FontSize',16);
    box off;
    xtickformat('$%,.0f')
    xlim([min_x,80])

    myAxes=findobj(s4,'Type','Axes');
    exportgraphics(myAxes,['National_Cot_per_case_' Xlabel_Add_TXT{:} '.png'],'Resolution',300);

    text(-0.15,1.05,'D','Fontsize',40,'Units','Normalized')

    print(gcf,['National_Measures_' Xlabel_Add_TXT{:} '.png'],'-dpng','-r300');
end


% print(gcf,['Preliminary_National.png'],'-dpng','-r300');
% 
% p1000=1-logncdf(1000,log_cases_baseline(1),log_cases_baseline(2));
% display(['Baseline: Probability outbreak over 1000 cases = ' num2str(p1000)]);
% 
% p1500=1-logncdf(1500,log_cases_baseline(1),log_cases_baseline(2));
% display(['Baseline: Probability outbreak over 1500 cases = ' num2str(p1500)]);
% 
% p2000=1-logncdf(2000,log_cases_baseline(1),log_cases_baseline(2));
% display(['Baseline: Probability outbreak over 2000 cases = ' num2str(p2000)]);
% 
% [~,Unvaccinated_Cases_County_Red_5,Vaccinated_Cases_County_Red_5,logn_Red_5]=Monte_Carlo_Incidence(0.05,NS);
% Unvaccinated_Cases_County_Red_5=squeeze(sum(Unvaccinated_Cases_County_Red_5,1));
% Vaccinated_Cases_County_Red_5=squeeze(sum(Vaccinated_Cases_County_Red_5,1));
% 
% [~,Unvaccinated_Cases_County_Red_10,Vaccinated_Cases_County_Red_10,logn_Red_10]=Monte_Carlo_Incidence(0.10,NS);
% Unvaccinated_Cases_County_Red_10=squeeze(sum(Unvaccinated_Cases_County_Red_10,1));
% Vaccinated_Cases_County_Red_10=squeeze(sum(Vaccinated_Cases_County_Red_10,1));
% 
% close all;
% 
% figure('units','normalized','outerposition',[0 0.2 1 0.6]);
% 
% Age_X=cell(18,1);
% subplot("Position",[0.075 0.175 0.875 0.8])
% for ii=1:18
%     if(ii<18)
%         Age_X{ii}=[num2str(5.*(ii-1)) char(8211) num2str(5.*(ii)-1)];
%     else
%         Age_X{ii}=[num2str(5.*(ii-1)) '+'];
%     end
%     temp_c=Unvaccinated_Cases_County_Baseline(ii,:);
%     temp_c(temp_c==0)=0.1;
%     log_temp=lognfit(temp_c);
%     log_temp=fmincon(@(x)logn_fit_zero(x,Unvaccinated_Cases_County_Baseline(ii,:)),log_temp);
% 
%     cc=linspace(0,logninv(0.995,log_temp(1),log_temp(2)),5001);
%     xx=lognpdf(cc,log_temp(1),log_temp(2));
%     xx=0.20.*xx./max(xx);
%     patch((ii-0.205)+[xx -flip(xx)],[cc flip(cc)],'r','LineStyle','none'); hold on
% 
%     temp_c=Vaccinated_Cases_County_Baseline(ii,:);
%     temp_c(temp_c==0)=0.1;
%     log_temp=lognfit(temp_c);
%     log_temp=fmincon(@(x)logn_fit_zero(x,Vaccinated_Cases_County_Baseline(ii,:)),log_temp);
% 
% 
% 
%     cc=linspace(0,logninv(0.995,log_temp(1),log_temp(2)),5001);
%     xx=lognpdf(cc,log_temp(1),log_temp(2));
%     xx=0.20.*xx./max(xx);
%     patch((ii+0.205)+[xx -flip(xx)],[cc flip(cc)],'b','LineStyle','none'); hold on
% end
% ylim([0 150])
% xlim([0.5 18.5]);
% set(gca,'LineWidth',2,'tickdir','out','FontSize',16,'XTick',[1:18],'XTickLabel',Age_X,'YTick',[0:25:150]);
% xlabel('Age group (years)','FontSize',18)
% ylabel('Annual measles cases','FontSize',18)
% legend({'Unvaccinated';'Vaccinated'},'FontSize',16)
% 
% print(gcf,['Preliminary_Age_Cases.png'],'-dpng','-r300');
% 
% figure('units','normalized','outerposition',[0 0.2 1 0.6]);
% 
% Age_X=cell(18,1);
% subplot("Position",[0.075 0.175 0.875 0.8])
% for ii=1:18
%     if(ii<18)
%         Age_X{ii}=[num2str(5.*(ii-1)) char(8211) num2str(5.*(ii)-1)];
%     else
%         Age_X{ii}=[num2str(5.*(ii-1)) '+'];
%     end
%     temp_c=Unvaccinated_Cases_County_Baseline(ii,:)+Vaccinated_Cases_County_Baseline(ii,:);
%     temp_c(temp_c==0)=0.1;
%     log_temp=lognfit(temp_c);
%     log_temp=fmincon(@(x)logn_fit_zero(x,Unvaccinated_Cases_County_Baseline(ii,:)+Vaccinated_Cases_County_Baseline(ii,:)),log_temp);
% 
%     cc=linspace(0,logninv(0.995,log_temp(1),log_temp(2)),5001);
%     xx=lognpdf(cc,log_temp(1),log_temp(2));
%     xx=0.125.*xx./max(xx);
%     patch((ii-0.3333)+[xx -flip(xx)],[cc flip(cc)],'k','LineStyle','none'); hold on
% 
%     temp_c=Unvaccinated_Cases_County_Red_5(ii,:)+Vaccinated_Cases_County_Red_5(ii,:);
%     temp_c(temp_c==0)=0.1;
%     log_temp=lognfit(temp_c);
%     log_temp=fmincon(@(x)logn_fit_zero(x,Unvaccinated_Cases_County_Red_5(ii,:)+Vaccinated_Cases_County_Red_5(ii,:)),log_temp);
% 
%     cc=linspace(0,logninv(0.995,log_temp(1),log_temp(2)),5001);
%     xx=lognpdf(cc,log_temp(1),log_temp(2));
%     xx=0.125.*xx./max(xx);
%     patch((ii)+[xx -flip(xx)],[cc flip(cc)],'b','LineStyle','none'); hold on
% 
%     temp_c=Unvaccinated_Cases_County_Red_10(ii,:)+Vaccinated_Cases_County_Red_10(ii,:);
%     temp_c(temp_c==0)=0.1;
%     log_temp=lognfit(temp_c);
%     log_temp=fmincon(@(x)logn_fit_zero(x,Unvaccinated_Cases_County_Red_10(ii,:)+Vaccinated_Cases_County_Red_10(ii,:)),log_temp);
% 
%     cc=linspace(0,logninv(0.995,log_temp(1),log_temp(2)),5001);
%     xx=lognpdf(cc,log_temp(1),log_temp(2));
%     xx=0.125.*xx./max(xx);
%     patch((ii+0.3333)+[xx -flip(xx)],[cc flip(cc)],'r','LineStyle','none'); hold on
% end
% ylim([0 150])
% xlim([0.5 18.5]);
% set(gca,'LineWidth',2,'tickdir','out','FontSize',16,'XTick',[1:18],'XTickLabel',Age_X,'YTick',[0:25:150]);
% xlabel('Age group (years)','FontSize',18)
% ylabel('Annual measles cases','FontSize',18)
% legend({'Baseline';'5% reduction'; '10% reduction'},'FontSize',16)
% 
% print(gcf,['Preliminary_Age_Cases_Reduction.png'],'-dpng','-r300');
% 
% 
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % County
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % C=({'#ffffff';
% %      '#fff5f0';
% % '#fee0d2';
% % '#fcbba1';
% % '#fc9272';
% % '#fb6a4a';
% % '#ef3b2c';
% % '#cb181d';
% % '#a50f15';
% % '#67000d'; ...
% % '#000000'});
% % 
% % 
% % C_Cases=hex2rgb(C);
% % x_avg_case=linspace(0,1,11);
% % 
% % 
% % states = shaperead('usastatelo', 'UseGeoCoords', true);
% % figure('units','normalized','outerposition',[0 0.075 1 1]);
% % 
% % subplot("Position",[0.9 0.05 0.01 0.9])
% % avg_c=linspace(0,1,1001);
% % CC_Cases=interp1(x_avg_case,C_Cases,avg_c);
% % 
% % for jj=1:1000
% %     patch([0 0 1 1],[avg_c(jj) avg_c(jj+1) avg_c(jj+1) avg_c(jj)],CC_Cases(jj,:),'LineStyle','None'); hold on;
% % end
% % 
% % plot([0 1],[0 0],'k','LineWidth',2)
% % plot([0 1],[1 1],'k','LineWidth',2)
% % plot([0 0],[0 1],'k','LineWidth',2)
% % plot([1 1],[0 1],'k','LineWidth',2)
% % 
% % text_v=[0:0.1:1];
% % for jj=1:length(text_v)
% %     if(jj==11)
% %         text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     else
% %         text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     end
% % end
% % text(5.45,0.5,'Average measles cases per year','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
% % axis off;
% % ylim([0 1])
% % 
% % 
% % S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
% % M=mean(Outbreak_County_VUR_10,2);
% % 
% % Avg_Mealses_County=zeros(length(S),1);
% % for ss=1:length(S)
% %     tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
% %     if(sum(tf)>0)
% %         Avg_Mealses_County(ss)=M(tf);
% %     end
% % end
% % 
% % Avg_Mealses_County(Avg_Mealses_County>1)=1;
% % 
% % 
% % NS=length(S);
% % 
% % 
% % 
% % 
% % 
% % CC_Cases=interp1(x_avg_case,C_Cases,Avg_Mealses_County);
% % CC_Cases(isnan(Avg_Mealses_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Mealses_County)),1);
% % CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});
% % 
% % ax1=usamap('conus');
% %  framem off; gridm off; mlabel off; plabel off;
% % ax1.Position=[-0.3,-2,0.6,0.6];
% % geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
% % geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 
% % 
% % ax1.Position=[-0.15,-0.15,1.2,1.2];
% % 
% % print(gcf,['Preliminary_County_Incidence_reduction_2_5.png'],'-dpng','-r300');
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % County
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % C=({'#ffffff';
% %      '#fff5f0';
% % '#fee0d2';
% % '#fcbba1';
% % '#fc9272';
% % '#fb6a4a';
% % '#ef3b2c';
% % '#cb181d';
% % '#a50f15';
% % '#67000d'; ...
% % '#000000'});
% % 
% % 
% % C_Cases=hex2rgb(C);
% % x_avg_case=linspace(0,1,11);
% % 
% % 
% % states = shaperead('usastatelo', 'UseGeoCoords', true);
% % figure('units','normalized','outerposition',[0 0.075 1 1]);
% % 
% % subplot("Position",[0.9 0.05 0.01 0.9])
% % avg_c=linspace(0,1,1001);
% % CC_Cases=interp1(x_avg_case,C_Cases,avg_c);
% % 
% % for jj=1:1000
% %     patch([0 0 1 1],[avg_c(jj) avg_c(jj+1) avg_c(jj+1) avg_c(jj)],CC_Cases(jj,:),'LineStyle','None'); hold on;
% % end
% % 
% % plot([0 1],[0 0],'k','LineWidth',2)
% % plot([0 1],[1 1],'k','LineWidth',2)
% % plot([0 0],[0 1],'k','LineWidth',2)
% % plot([1 1],[0 1],'k','LineWidth',2)
% % 
% % text_v=[0:0.1:1];
% % for jj=1:length(text_v)
% %     if(jj==11)
% %         text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     else
% %         text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     end
% % end
% % text(5.45,0.5,'Average measles cases per year','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
% % axis off;
% % ylim([0 1])
% % 
% % 
% % S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
% % M=mean(Outbreak_County,2);
% % 
% % Avg_Mealses_County=zeros(length(S),1);
% % for ss=1:length(S)
% %     tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
% %     if(sum(tf)>0)
% %         Avg_Mealses_County(ss)=M(tf);
% %     end
% % end
% % 
% % Avg_Mealses_County(Avg_Mealses_County>1)=1;
% % 
% % 
% % NS=length(S);
% % 
% % 
% % 
% % 
% % 
% % CC_Cases=interp1(x_avg_case,C_Cases,Avg_Mealses_County);
% % CC_Cases(isnan(Avg_Mealses_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Mealses_County)),1);
% % CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});
% % 
% % ax1=usamap('conus');
% %  framem off; gridm off; mlabel off; plabel off;
% % ax1.Position=[-0.3,-2,0.6,0.6];
% % geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
% % geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 
% % 
% % ax1.Position=[-0.15,-0.15,1.2,1.2];
% % 
% % print(gcf,['Preliminary_County_Incidence_Baseline.png'],'-dpng','-r300');
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % County
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % C=({'#ffffff';
% %      '#f7f4f9';
% % '#e7e1ef';
% % '#d4b9da';
% % '#c994c7';
% % '#df65b0';
% % '#e7298a';
% % '#ce1256';
% % '#980043';
% % '#67001f'; ...
% % '#000000'});
% % 
% % 
% % C_Risk=hex2rgb(C);
% % x_risk=linspace(0,0.1,11);
% % 
% % 
% % states = shaperead('usastatelo', 'UseGeoCoords', true);
% % figure('units','normalized','outerposition',[0 0.075 1 1]);
% % 
% % subplot("Position",[0.9 0.05 0.01 0.9])
% % avg_risk=linspace(0,0.1,1001);
% % CC_Risk=interp1(x_risk,C_Risk,avg_risk);
% % 
% % for jj=1:1000
% %     patch([0 0 1 1],[avg_risk(jj) avg_risk(jj+1) avg_risk(jj+1) avg_risk(jj)],CC_Risk(jj,:),'LineStyle','None'); hold on;
% % end
% % 
% % plot([0 1],[0 0],'k','LineWidth',2)
% % plot([0 1],0.1.*[1 1],'k','LineWidth',2)
% % plot([0 0],0.1.*[0 1],'k','LineWidth',2)
% % plot([1 1],0.1.*[0 1],'k','LineWidth',2)
% % 
% % text_v=[0:0.01:0.1];
% % for jj=1:length(text_v)
% %     if(jj==11)
% %         text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     else
% %         text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     end
% % end
% % text(5.45,0.05,'Probability of at least one measles case','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
% % axis off;
% % ylim([0 0.1])
% % 
% % 
% % S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
% % M=1-p_zero_10;
% % 
% % Avg_Risk_County=zeros(length(S),1);
% % for ss=1:length(S)
% %     tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
% %     if(sum(tf)>0)
% %         Avg_Risk_County(ss)=M(tf);
% %     end
% % end
% % Avg_Risk_County(Avg_Risk_County>0.1)=0.1;
% % NS=length(S);
% % 
% % 
% % 
% % 
% % 
% % CC_Cases=interp1(x_risk,C_Risk,Avg_Risk_County);
% % CC_Cases(isnan(Avg_Risk_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Risk_County)),1);
% % CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});
% % 
% % ax1=usamap('conus');
% %  framem off; gridm off; mlabel off; plabel off;
% % ax1.Position=[-0.3,-2,0.6,0.6];
% % geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
% % geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 
% % 
% % ax1.Position=[-0.15,-0.15,1.2,1.2];
% % 
% % print(gcf,['Preliminary_County_Risk_reduction_2_5.png'],'-dpng','-r300');
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% % % County
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % C=({'#ffffff';
% %      '#f7f4f9';
% % '#e7e1ef';
% % '#d4b9da';
% % '#c994c7';
% % '#df65b0';
% % '#e7298a';
% % '#ce1256';
% % '#980043';
% % '#67001f'; ...
% % '#000000'});
% % 
% % 
% % C_Risk=hex2rgb(C);
% % x_risk=linspace(0,0.1,11);
% % 
% % 
% % states = shaperead('usastatelo', 'UseGeoCoords', true);
% % figure('units','normalized','outerposition',[0 0.075 1 1]);
% % 
% % subplot("Position",[0.9 0.05 0.01 0.9])
% % avg_risk=linspace(0,0.1,1001);
% % CC_Risk=interp1(x_risk,C_Risk,avg_risk);
% % 
% % for jj=1:1000
% %     patch([0 0 1 1],[avg_risk(jj) avg_risk(jj+1) avg_risk(jj+1) avg_risk(jj)],CC_Risk(jj,:),'LineStyle','None'); hold on;
% % end
% % 
% % plot([0 1],[0 0],'k','LineWidth',2)
% % plot([0 1],0.1.*[1 1],'k','LineWidth',2)
% % plot([0 0],0.1.*[0 1],'k','LineWidth',2)
% % plot([1 1],0.1.*[0 1],'k','LineWidth',2)
% % 
% % text_v=[0:0.01:0.1];
% % for jj=1:length(text_v)
% %     if(jj==11)
% %         text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     else
% %         text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
% %     end
% % end
% % text(5.45,0.05,'Probability of at least one measles case','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
% % axis off;
% % ylim([0 0.1])
% % 
% % 
% % S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
% % M=1-p_zero;
% % 
% % Avg_Risk_County=zeros(length(S),1);
% % for ss=1:length(S)
% %     tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
% %     if(sum(tf)>0)
% %         Avg_Risk_County(ss)=M(tf);
% %     end
% % end
% % Avg_Risk_County(Avg_Risk_County>0.1)=0.1;
% % NS=length(S);
% % 
% % 
% % 
% % 
% % 
% % CC_Cases=interp1(x_risk,C_Risk,Avg_Risk_County);
% % CC_Cases(isnan(Avg_Risk_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Risk_County)),1);
% % CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});
% % 
% % ax1=usamap('conus');
% %  framem off; gridm off; mlabel off; plabel off;
% % ax1.Position=[-0.3,-2,0.6,0.6];
% % geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
% % geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 
% % 
% % ax1.Position=[-0.15,-0.15,1.2,1.2];
% % 
% print(gcf,['Preliminary_County_Risk_Baseline.png'],'-dpng','-r300');