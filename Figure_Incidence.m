function Figure_Incidence(Scenario_Plot)
    close all;
    
    load('Turncated_Negative_Binomial_Parameter.mat');
    F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
    
    Age_Reduction=[true(1,4) false(1,1)];
    NS=2.5.*10^3;
    National_Reduction=0;
    [p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();
    [Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization]=Measles_Outbreak_Cost();
    
    [Outbreak_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Incidence(F_NB,National_Reduction,Age_Reduction,NS,Scenario_Plot);
    
    Unvaccinated_Cases_County_Baseline=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
    Vaccinated_Cases_County_Baseline=squeeze(sum(Vaccinated_Cases_County_Baseline,1));
    
    Hospitalizations_Baseline=p_H_Unvaccinated*Unvaccinated_Cases_County_Baseline+p_H_Vaccinated*Vaccinated_Cases_County_Baseline;
    
    Total_Contacts_Baseline=squeeze(sum(Total_Contacts_Baseline,[1 2]));
    Cost_Vac_Age=[Cost_per_Vaccine_dose_VFC.*ones(1,4) Cost_per_Vaccine_dose_Private.*ones(1,14)];
    Cost_Vaccination_Contacts=Cost_Vac_Age*squeeze(sum(Unvaccinated_Contacts_Baseline,1));
    Cost_Case_Medical=Cost_per_Hospitalization.*Hospitalizations_Baseline+Cost_per_Non_Hospitalization.*(sum(Outbreak_Cases_County,1)-Hospitalizations_Baseline);
    Cost_Baseline=Cost_per_Contact.*Total_Contacts_Baseline+Cost_Vaccination_Contacts(:)+Cost_Case_Medical(:);
    
    
    temp_c=sum(Outbreak_Cases_County,1);
    pd=fitdist(temp_c(:),'Kernel','Support','positive');
    % temp_c(temp_c==0)=1-mean(temp_c==0);
    % log_temp=lognfit(temp_c);
    % log_cases_baseline=fmincon(@(x)logn_fit_zero(x,sum(Outbreak_Cases_County,1)),log_temp);
    
    temp_c=Hospitalizations_Baseline;
    temp_c(temp_c==0)=1-mean(temp_c==0);
    log_temp=lognfit(temp_c);
    log_hospital_baseline=fmincon(@(x)logn_fit_zero(x,Hospitalizations_Baseline),log_temp);
    
    temp_c=Cost_Baseline;
    temp_c(temp_c==0)=1-mean(temp_c==0);
    log_temp=lognfit(temp_c);
    log_cost_baseline=fmincon(@(x)logn_fit_zero(x,Cost_Baseline),log_temp);
    
    
    Cost_per_Case=Cost_Baseline./(sum(Outbreak_Cases_County,1)');
    temp_c=Cost_per_Case;
    temp_c(temp_c==0)=1-mean(temp_c==0);
    log_temp=lognfit(temp_c);
    log_cost_per_case_baseline=fmincon(@(x)logn_fit_zero(x,Cost_per_Case),log_temp);

    
    figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    
    subplot("Position",[0.075 0.585 0.4 0.375])
    cases=linspace(0,5000,10^4);
    plot(cases,pdf(pd,cases),'k');%,log_cases_baseline(1),log_cases_baseline(2)),'k');%,cases,lognpdf(cases,log_cases_5(1),log_cases_5(2)),'b',cases,lognpdf(cases,log_cases_10(1),log_cases_10(2)),'r','LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[0:500:5000],'Fontsize',16)
    xlabel('Annual measles cases','FontSize',18)
    ylabel('Probability density','FontSize',18)
    %legend({'Baseline','5% reduction','10% reduction'},'FontSize',16);
    box off;
    subplot("Position",[0.575 0.585 0.4 0.375])
    hospital=linspace(0,1750,10^4);
    plot(hospital,lognpdf(hospital,log_hospital_baseline(1),log_hospital_baseline(2)),'k');%,hospital,lognpdf(hospital,log_hospital_5(1),log_hospital_5(2)),'b',hospital,lognpdf(hospital,log_hospital_10(1),log_hospital_10(2)),'r','LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',[0:250:1750],'Fontsize',16)
    xlabel('Annual measles hospitalizations','FontSize',18)
    ylabel('Probability density','FontSize',18)
    %legend({'Baseline','5% reduction','10% reduction'},'FontSize',16);
    box off;
    
    subplot("Position",[0.075 0.105 0.4 0.375])
    cost=10.^7.*linspace(0,10,10^4);
    plot(cost,lognpdf(cost,log_cost_baseline(1),log_cost_baseline(2)),'k');%,cost,lognpdf(cost,log_cost_5(1),log_cost_5(2)),'b',cost,lognpdf(cost,log_cost_10(1),log_cost_10(2)),'r','LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',10.^7.*[0:1:10],'Fontsize',16)
    xlabel('Total costs','FontSize',18)
    ylabel('Probability density','FontSize',18)
    %legend({'Baseline','5% reduction','10% reduction'},'FontSize',16);
    box off;
    subplot("Position",[0.575 0.105 0.4 0.375])
    cost_per_case=10.^4.*linspace(2.3,3,10^4);
    plot(cost_per_case,lognpdf(cost_per_case,log_cost_per_case_baseline(1),log_cost_per_case_baseline(2)),'k');%,cost_per_case,lognpdf(cost_per_case,log_cost_per_case_5(1),log_cost_per_case_5(2)),'b',cost_per_case,lognpdf(cost_per_case,log_cost_per_case_10(1),log_cost_per_case_10(2)),'r','LineWidth',2);
    set(gca,'LineWidth',2,'Tickdir','out','XTick',10.^4.*[2.3:0.1:3],'Fontsize',16)
    xlabel('Cost per case','FontSize',18)
    ylabel('Probability density','FontSize',18)
    %legend({'Baseline','5% reduction','10% reduction'},'FontSize',16);
    box off;
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
% % print(gcf,['Preliminary_County_Incidence_Reduction_10.png'],'-dpng','-r300');
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
% % print(gcf,['Preliminary_County_Risk_Reduction_10.png'],'-dpng','-r300');
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