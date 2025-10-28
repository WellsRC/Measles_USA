% function Figure_Comparison_County(National_Reduction,Scenario,Age_0_to_6)
% 
National_Reduction=0.01;
Scenario='Sample_2025';
Age_0_to_6=true;

close all;

if(~Age_0_to_6)
    Age_Reduction=[true(1,1) false(1,4)];
else
    Age_Reduction=[true(1,2) false(1,3)];
end
N_Samp=2.5.*10^3;

[cases_baseline,hospital_baseline,cost_baseline,cost_per_case_baseline]=County_Outcome_Central_Measure(0,Age_Reduction,N_Samp,Scenario,Age_0_to_6);
[cases_reduction,hospital_reduction,cost_reduction,cost_per_case_reduction]=County_Outcome_Central_Measure(National_Reduction,Age_Reduction,N_Samp,Scenario,Age_0_to_6);
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);

FN_Age_Class={'Age_0_to_4','_Age_5_to_9','_Age_10_to_14','_Age_15_to_19','_Age_20_to_24'};
FN_Age_Class=FN_Age_Class(Age_Reduction);

if(Age_0_to_6)
    load(['National_Reduction=' num2str(100*0) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction')
else
    load(['National_Reduction=' num2str(100*0) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')
end

V_Baseline=NaN.*zeros(length(S),1);
V_Reduction=NaN.*zeros(length(S),1);

Case_Baseline=NaN.*zeros(length(S),1);
Case_Reduction=NaN.*zeros(length(S),1);

Hospital_Baseline=NaN.*zeros(length(S),1);
Hospital_Reduction=NaN.*zeros(length(S),1);

Cost_Baseline=NaN.*zeros(length(S),1);
Cost_Reduction=NaN.*zeros(length(S),1);

Cost_per_case_Baseline=NaN.*zeros(length(S),1);
Cost_per_case_Reduction=NaN.*zeros(length(S),1);

for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        if(~Age_0_to_6)
            V_Baseline(cc,:)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
        else
            V_Baseline(cc,:)=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,1)).*table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,2)).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)))./(table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)));
        end
        Case_Baseline(cc,:)=10.^4.*cases_baseline(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));
        Case_Reduction(cc,:)=10.^4.*cases_reduction(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));

        Hospital_Baseline(cc,:)=10.^6.*hospital_baseline(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));
        Hospital_Reduction(cc,:)=10.^6.*hospital_reduction(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));

        Cost_Baseline(cc,:)=10.^4.*cost_baseline(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));
        Cost_Reduction(cc,:)=10.^4.*cost_reduction(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));

        Cost_per_case_Baseline(cc,:)=cost_per_case_baseline(tf);
        Cost_per_case_Reduction(cc,:)=cost_per_case_reduction(tf);
    end
end

if(Age_0_to_6)
    load(['National_Reduction=' num2str(100*National_Reduction) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction')
else
    load(['National_Reduction=' num2str(100*National_Reduction) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')
end

for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
     if(sum(tf)>0)
        if(~Age_0_to_6)
            V_Reduction(cc,:)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
        else
            V_Reduction(cc,:)=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,1)).*table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,2)).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)))./(table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)));
        end
    end
end

V_Diff=V_Baseline-V_Reduction;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Plot Vaccine Uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 C=flip({'#ffffff';
     '#ffffcc';
'#ffeda0';
'#fed976';
'#feb24c';
'#fd8d3c';
'#fc4e2a';
'#e31a1c';
'#bd0026';
'#800026'; ...
'#000000'});


C_vac=hex2rgb(C);
x_vac=linspace(0.05,1,11);



 C=flip({'#b2182b';
     '#ffffff';
     '#2166ac';});


C_diff=hex2rgb(C);
x_diff=[0 National_Reduction, ceil(100.*(max(V_Diff(:))))./100];

states = shaperead('usastatelo', 'UseGeoCoords', true);

figure('units','normalized','outerposition',[0.15 0.075 0.7 0.4]);

subplot("Position",[0.0225 0.135 0.28 0.025])
dv=linspace(0.05,1,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch(100.*[dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([5 5],[0 1],'k','LineWidth',1)
plot([100 100],[0 1],'k','LineWidth',1)
plot([5 100],[0 0],'k','LineWidth',1)
plot([5 100],[1 1],'k','LineWidth',1)

text_v=[5 20:20:100];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj)) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj)) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(52.5,-3.9,'MMR vaccine uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
text(0.405,33.5,'A','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle')
text(52.5,33.5,'Baseline','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
xlim([5 100])

subplot("Position",[0.36 0.135 0.28 0.025])
dv=linspace(0.05,1,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch(100.*[dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([5 5],[0 1],'k','LineWidth',1)
plot([100 100],[0 1],'k','LineWidth',1)
plot([5 100],[0 0],'k','LineWidth',1)
plot([5 100],[1 1],'k','LineWidth',1)

text_v=[5 20:20:100];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj)) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj)) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(52.5,-3.9,'MMR vaccine uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
text(0.405,33.5,'B','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle')
text(52.5,33.5,'Reduction','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
xlim([5 100])


subplot("Position",[0.7 0.135 0.28 0.025])
dv=linspace(0.0,max(x_diff),1001);
CC_Diff=interp1(x_diff,C_diff,dv);

for jj=1:1000
    patch(100.*[dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Diff(jj,:),'LineStyle','None'); hold on;
end

plot([0 0],[0 1],'k','LineWidth',1)
plot(max(x_diff).*[100 100],[0 1],'k','LineWidth',1)
plot(max(x_diff).*[0 100],[0 0],'k','LineWidth',1)
plot(max(x_diff).*[0 100],[1 1],'k','LineWidth',1)

text_v=round(linspace(0,100.*max(x_diff),6));
for jj=1:length(text_v)
    text(text_v(jj),-1.2,[num2str(text_v(jj)) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
end
text(50.*max(x_diff),-3.9,'Absolute reduction in uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
text(-0.05,33.5,'C','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Difference','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
xlim([0 100.*max(x_diff)])

V_Baseline(V_Baseline<0.05)=0.05;
V_Reduction(V_Reduction<0.05)=0.05;

NS=length(S);

CC_Vac=interp1(x_vac,C_vac,V_Baseline(:));
CC_Vac(isnan(V_Baseline(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Baseline(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,V_Reduction(:));
CC_Vac(isnan(V_Reduction(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Reduction(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax2=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax2.Position=[1.7,0.4,0.6,0.6];
geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,V_Diff(:));
CC_Diff(isnan(V_Diff(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Diff(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax3=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax3.Position=[-0.3,-0.1,0.6,0.6];

geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); 


ax1.Position=[-0.32, -0.015 ,1,1];
ax2.Position=[0.02, -0.015 ,1,1];
ax3.Position=[0.35, -0.015 ,1,1];

print(gcf,['Figure_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 C={'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d';
'#000000'};


C_vac=hex2rgb(C);

mx=0.01;
MX=round(100.*max(prctile(Case_Baseline,[99]),prctile(Case_Reduction,[99])))./100;
x_vac=linspace(mx,MX,11);


dC=Case_Reduction-Case_Baseline;
 C=({'#ffffff';
'#fde0dd';
'#fcc5c0';
'#fa9fb5';
'#f768a1';
'#dd3497';
'#ae017e';
'#7a0177';
'#49006a';
'#000000'});


C_diff=hex2rgb(C);
x_diff=linspace(0,ceil(100.*prctile(dC,[99]))./100,size(C_diff,1));

states = shaperead('usastatelo', 'UseGeoCoords', true);

figure('units','normalized','outerposition',[0.15 0.075 0.7 0.4]);

subplot("Position",[0.0225 0.135 0.28 0.025])
dv=linspace(mx,MX,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([mx mx],[0 1],'k','LineWidth',1)
plot([MX MX],[0 1],'k','LineWidth',1)
plot([mx MX],[0 0],'k','LineWidth',1)
plot([mx MX],[1 1],'k','LineWidth',1)

text_v=[mx 1:(floor(MX)-1) MX];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    elseif(jj==length(text_v))
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj)) ],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Cases per 10,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(-0.05,33.5,'A','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Baseline','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
xlim([mx MX])

subplot("Position",[0.36 0.135 0.28 0.025])
dv=linspace(mx,MX,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([mx mx],[0 1],'k','LineWidth',1)
plot([MX MX],[0 1],'k','LineWidth',1)
plot([mx MX],[0 0],'k','LineWidth',1)
plot([mx MX],[1 1],'k','LineWidth',1)

text_v=[mx 1:(floor(MX)-1) MX];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    elseif(jj==length(text_v))
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Cases per 10,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(-0.05,33.5,'B','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Reduction','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
xlim([mx MX])


subplot("Position",[0.7 0.135 0.28 0.025])
dv=linspace(0.0,max(x_diff),1001);
CC_Diff=interp1(x_diff,C_diff,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Diff(jj,:),'LineStyle','None'); hold on;
end

plot([0 0],[0 1],'k','LineWidth',1)
plot(max(x_diff).*[1 1],[0 1],'k','LineWidth',1)
plot(max(x_diff).*[0 1],[0 0],'k','LineWidth',1)
plot(max(x_diff).*[0 1],[1 1],'k','LineWidth',1)

text_v=round(linspace(0,max(x_diff),6),2);
for jj=1:length(text_v)
    if(jj<length(text_v))
        text(text_v(jj),-1.2,[num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Additional cases per 10,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
text(-0.05,33.5,'C','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Difference','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
xlim([0 max(x_diff)])

Case_Baseline(Case_Baseline<mx)=mx;
Case_Reduction(Case_Reduction<mx)=mx;
Case_Baseline(Case_Baseline>MX)=MX;
Case_Reduction(Case_Reduction>MX)=MX;
dC(dC>max(x_diff))=max(x_diff);

NS=length(S);

CC_Vac=interp1(x_vac,C_vac,Case_Baseline(:));
CC_Vac(isnan(Case_Baseline(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(Case_Baseline(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,Case_Reduction(:));
CC_Vac(isnan(Case_Reduction(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(Case_Reduction(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax2=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax2.Position=[1.7,0.4,0.6,0.6];
geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,dC(:));
CC_Diff(isnan(dC(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(dC(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax3=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax3.Position=[-0.3,-0.1,0.6,0.6];

geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); 


ax1.Position=[-0.32, -0.015 ,1,1];
ax2.Position=[0.02, -0.015 ,1,1];
ax3.Position=[0.35, -0.015 ,1,1];

print(gcf,['Figure_Cases_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hospitalizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 C={'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d';
'#000000'};


C_vac=hex2rgb(C);

mx=round(100.*min(prctile(Hospital_Baseline,[1]),prctile(Hospital_Reduction,[1])))./100;
MX=round(100.*max(prctile(Hospital_Baseline,[99]),prctile(Hospital_Reduction,[99])))./100;
x_vac=linspace(mx,MX,11);


dC=Hospital_Reduction-Hospital_Baseline;
 C=({'#ffffff';
'#fde0dd';
'#fcc5c0';
'#fa9fb5';
'#f768a1';
'#dd3497';
'#ae017e';
'#7a0177';
'#49006a';
'#000000'});


C_diff=hex2rgb(C);
x_diff=linspace(floor(100.*prctile(dC,[1]))./100,ceil(100.*prctile(dC,[99]))./100,size(C_diff,1));

states = shaperead('usastatelo', 'UseGeoCoords', true);

figure('units','normalized','outerposition',[0.15 0.075 0.7 0.4]);

subplot("Position",[0.0225 0.135 0.28 0.025])
dv=linspace(mx,MX,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([mx mx],[0 1],'k','LineWidth',1)
plot([MX MX],[0 1],'k','LineWidth',1)
plot([mx MX],[0 0],'k','LineWidth',1)
plot([mx MX],[1 1],'k','LineWidth',1)

text_v=[mx 40:40:MX];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    elseif(jj==length(text_v))
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj)) ],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Hospitalizations per 1,000,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(-0.05,33.5,'A','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Baseline','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
xlim([mx MX])

subplot("Position",[0.36 0.135 0.28 0.025])
dv=linspace(mx,MX,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([mx mx],[0 1],'k','LineWidth',1)
plot([MX MX],[0 1],'k','LineWidth',1)
plot([mx MX],[0 0],'k','LineWidth',1)
plot([mx MX],[1 1],'k','LineWidth',1)

text_v=[mx 40:40:MX];
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    elseif(jj==length(text_v))
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,[num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Hospitalizations per 1,000,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(-0.05,33.5,'B','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Reduction','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
xlim([mx MX])


subplot("Position",[0.7 0.135 0.28 0.025])
dv=linspace(min(x_diff),max(x_diff),1001);
CC_Diff=interp1(x_diff,C_diff,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Diff(jj,:),'LineStyle','None'); hold on;
end

plot([min(x_diff) min(x_diff)],[0 1],'k','LineWidth',1)
plot([max(x_diff) max(x_diff)],[0 1],'k','LineWidth',1)
plot([min(x_diff) max(x_diff)],[0 0],'k','LineWidth',1)
plot([min(x_diff) max(x_diff)],[1 1],'k','LineWidth',1)

text_v=round(linspace(min(x_diff),max(x_diff),6),2);
for jj=1:length(text_v)
    if(jj==1)
        text(text_v(jj),-1.2,['\leq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    elseif(jj<length(text_v))
        text(text_v(jj),-1.2,[num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    else
        text(text_v(jj),-1.2,['\geq ' num2str(text_v(jj))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
    end
end
text(0.5,-3.9,'Additional hospitalizations per 1,000,000','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
axis off;
text(-0.05,33.5,'C','FontSize',28,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
text(0.5,33.5,'Difference','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized')
xlim([min(x_diff) max(x_diff)])

Hospital_Baseline(Hospital_Baseline<mx)=mx;
Hospital_Reduction(Hospital_Reduction<mx)=mx;
Hospital_Baseline(Hospital_Baseline>MX)=MX;
Hospital_Reduction(Hospital_Reduction>MX)=MX;
dC(dC>max(x_diff))=max(x_diff);
dC(dC<min(x_diff))=min(x_diff);

NS=length(S);

CC_Vac=interp1(x_vac,C_vac,Hospital_Baseline(:));
CC_Vac(isnan(Hospital_Baseline(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(Hospital_Baseline(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,Hospital_Reduction(:));
CC_Vac(isnan(Hospital_Reduction(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(Hospital_Reduction(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax2=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax2.Position=[1.7,0.4,0.6,0.6];
geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,dC(:));
CC_Diff(isnan(dC(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(dC(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax3=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax3.Position=[-0.3,-0.1,0.6,0.6];

geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); 


ax1.Position=[-0.32, -0.015 ,1,1];
ax2.Position=[0.02, -0.015 ,1,1];
ax3.Position=[0.35, -0.015 ,1,1];

print(gcf,['Figure_Hospitalizations_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');

