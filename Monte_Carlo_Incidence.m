clear;
clc;

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')
load('Baseline_Estimate_Measles_Incidence.mat',"beta_seed","k_nbin","beta_j","lambda_g",'k_mealses');


Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');
Imported_Case=zeros(length(County_Data.County),1);
for cc=1:length(Imported_Case)
    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'imported') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Imported_Case(cc)=Measles_Cases.case_count(t_f);
    end
end



load('County_Matrix_Gravity_Covariate.mat',"Gravity_Model",'County_GEOID');

indx_G=zeros(length(County_Data.County),1);
for cc=1:length(County_Data.County)
   indx_G(cc)=find(strcmp(County_GEOID,County_Data.GEOID{cc}));
end

G=Gravity_Model(indx_G,:);
G=G(:,indx_G);


Reff_Seed=interp1(County_Data.beta_j',County_Data.R_eff',beta_seed)';
Reff=interp1(County_Data.beta_j',County_Data.R_eff',beta_j)';

Case_Count=interp1(County_Data.beta_j',County_Data.Final_Size_Est',beta_j)';

Case_Count(Reff<1)=1./(1-Reff(Reff<1));

GM=G*Case_Count(:)./(length(Reff(:))-1);

p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^(Imported_Case+exp(lambda_g.*GM));

p_zero=zeros(size(p_c));

for ss=1:length(Reff)
    if(Reff(ss)>1)
        p_zero(ss)=p_c(ss)+(1-p_c(ss)).*nbinpdf(0,k_nbin,k_nbin./(k_nbin+Case_Count(ss)));
    else
        p_zero(ss)=p_c(ss)+(1-p_c(ss)).*Chain_Size_Distribution(1,Reff(ss),k_mealses);
    end
end

[Outbreak_County]=Monte_Carlo_Outbreak_County(p_c,Case_Count,k_nbin,Reff,k_mealses,10^4);
NOB=sum(Outbreak_County,1)';
NOB(NOB==0)=10^(-16);

logn_baseline=lognfit(NOB);

load(['National_Reduction_5.mat'],'County_Data_Vaccine_Reduction');
Reff_Seed=County_Data_Vaccine_Reduction.R_eff_Seed;
Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

Case_Count(Reff<1)=1./(1-Reff(Reff<1));

GM=G*Case_Count(:)./(length(Reff(:))-1);

k_mealses=0.23; 
p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^(Imported_Case+exp(lambda_g.*GM));


p_zero_5=zeros(size(p_c));

for ss=1:length(Reff)
    if(Reff(ss)>1)
        p_zero_5(ss)=p_c(ss)+(1-p_c(ss)).*nbinpdf(0,k_nbin,k_nbin./(k_nbin+Case_Count(ss)));
    else
        p_zero_5(ss)=p_c(ss)+(1-p_c(ss)).*Chain_Size_Distribution(1,Reff(ss),k_mealses);
    end
end


[Outbreak_County_VUR_5]=Monte_Carlo_Outbreak_County(p_c,Case_Count,k_nbin,Reff,k_mealses,10^4);
NOB=sum(Outbreak_County_VUR_5,1)';
NOB(NOB==0)=10^(-16);

logn_reduced=lognfit(NOB);


load(['National_Reduction_10.mat'],'County_Data_Vaccine_Reduction');
Reff_Seed=County_Data_Vaccine_Reduction.R_eff_Seed;
Reff=County_Data_Vaccine_Reduction.R_eff;
Case_Count=County_Data_Vaccine_Reduction.Final_Size_Est;

Case_Count(Reff<1)=1./(1-Reff(Reff<1));

GM=G*Case_Count(:)./(length(Reff(:))-1);

k_mealses=0.23; 
p_c=nbinpdf(0,k_mealses,k_mealses./(k_mealses+Reff_Seed)).^(Imported_Case+exp(lambda_g.*GM));

p_zero_10=zeros(size(p_c));

for ss=1:length(Reff)
    if(Reff(ss)>1)
        p_zero_10(ss)=p_c(ss)+(1-p_c(ss)).*nbinpdf(0,k_nbin,k_nbin./(k_nbin+Case_Count(ss)));
    else
        p_zero_10(ss)=p_c(ss)+(1-p_c(ss)).*Chain_Size_Distribution(1,Reff(ss),k_mealses);
    end
end

[Outbreak_County_VUR_10]=Monte_Carlo_Outbreak_County(p_c,Case_Count,k_nbin,Reff,k_mealses,10^4);
NOB=sum(Outbreak_County_VUR_10,1)';
NOB(NOB==0)=10^(-16);

logn_reduced_10=lognfit(NOB);

x=linspace(0,2500,10001);
close all;

figure(1)
plot(x,lognpdf(x,logn_baseline(1),logn_baseline(2)),'k',x,lognpdf(x,logn_reduced(1),logn_reduced(2)),'b-.',x,lognpdf(x,logn_reduced_10(1),logn_reduced_10(2)),'r-','LineWidth',2);
set(gca,'LineWidth',2,'TickDir','out','XTick',[0:250:2500],'FontSize',14);
xlim([0 2500])
xlabel('Annual measles incidence','FontSize',18)
ylabel('Probability density','FontSize',18)
legend({'Baseline';'5% reduction (0-9 year olds)';'10% reduction (0-9 year olds)'},'Fontsize',14)
box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C=({'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d'; ...
'#000000'});


C_Cases=hex2rgb(C);
x_avg_case=linspace(0,1,11);


states = shaperead('usastatelo', 'UseGeoCoords', true);
figure('units','normalized','outerposition',[0 0.075 1 1]);

subplot("Position",[0.9 0.05 0.01 0.9])
avg_c=linspace(0,1,1001);
CC_Cases=interp1(x_avg_case,C_Cases,avg_c);

for jj=1:1000
    patch([0 0 1 1],[avg_c(jj) avg_c(jj+1) avg_c(jj+1) avg_c(jj)],CC_Cases(jj,:),'LineStyle','None'); hold on;
end

plot([0 1],[0 0],'k','LineWidth',2)
plot([0 1],[1 1],'k','LineWidth',2)
plot([0 0],[0 1],'k','LineWidth',2)
plot([1 1],[0 1],'k','LineWidth',2)

text_v=[0:0.1:1];
for jj=1:length(text_v)
    if(jj==11)
        text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    else
        text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    end
end
text(5.45,0.5,'Average measles cases per year','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
ylim([0 1])


S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
M=mean(Outbreak_County_VUR_10,2);

Avg_Mealses_County=zeros(length(S),1);
for ss=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
    if(sum(tf)>0)
        Avg_Mealses_County(ss)=M(tf);
    end
end

Avg_Mealses_County(Avg_Mealses_County>1)=1;


NS=length(S);





CC_Cases=interp1(x_avg_case,C_Cases,Avg_Mealses_County);
CC_Cases(isnan(Avg_Mealses_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Mealses_County)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.15,-0.15,1.2,1.2];

print(gcf,['Preliminary_County_Incidence_Reduction_10.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C=({'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d'; ...
'#000000'});


C_Cases=hex2rgb(C);
x_avg_case=linspace(0,1,11);


states = shaperead('usastatelo', 'UseGeoCoords', true);
figure('units','normalized','outerposition',[0 0.075 1 1]);

subplot("Position",[0.9 0.05 0.01 0.9])
avg_c=linspace(0,1,1001);
CC_Cases=interp1(x_avg_case,C_Cases,avg_c);

for jj=1:1000
    patch([0 0 1 1],[avg_c(jj) avg_c(jj+1) avg_c(jj+1) avg_c(jj)],CC_Cases(jj,:),'LineStyle','None'); hold on;
end

plot([0 1],[0 0],'k','LineWidth',2)
plot([0 1],[1 1],'k','LineWidth',2)
plot([0 0],[0 1],'k','LineWidth',2)
plot([1 1],[0 1],'k','LineWidth',2)

text_v=[0:0.1:1];
for jj=1:length(text_v)
    if(jj==11)
        text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    else
        text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    end
end
text(5.45,0.5,'Average measles cases per year','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
ylim([0 1])


S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
M=mean(Outbreak_County,2);

Avg_Mealses_County=zeros(length(S),1);
for ss=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
    if(sum(tf)>0)
        Avg_Mealses_County(ss)=M(tf);
    end
end

Avg_Mealses_County(Avg_Mealses_County>1)=1;


NS=length(S);





CC_Cases=interp1(x_avg_case,C_Cases,Avg_Mealses_County);
CC_Cases(isnan(Avg_Mealses_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Mealses_County)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.15,-0.15,1.2,1.2];

print(gcf,['Preliminary_County_Incidence_Baseline.png'],'-dpng','-r300');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C=({'#ffffff';
     '#f7f4f9';
'#e7e1ef';
'#d4b9da';
'#c994c7';
'#df65b0';
'#e7298a';
'#ce1256';
'#980043';
'#67001f'; ...
'#000000'});


C_Risk=hex2rgb(C);
x_risk=linspace(0,0.1,11);


states = shaperead('usastatelo', 'UseGeoCoords', true);
figure('units','normalized','outerposition',[0 0.075 1 1]);

subplot("Position",[0.9 0.05 0.01 0.9])
avg_risk=linspace(0,0.1,1001);
CC_Risk=interp1(x_risk,C_Risk,avg_risk);

for jj=1:1000
    patch([0 0 1 1],[avg_risk(jj) avg_risk(jj+1) avg_risk(jj+1) avg_risk(jj)],CC_Risk(jj,:),'LineStyle','None'); hold on;
end

plot([0 1],[0 0],'k','LineWidth',2)
plot([0 1],0.1.*[1 1],'k','LineWidth',2)
plot([0 0],0.1.*[0 1],'k','LineWidth',2)
plot([1 1],0.1.*[0 1],'k','LineWidth',2)

text_v=[0:0.01:0.1];
for jj=1:length(text_v)
    if(jj==11)
        text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    else
        text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    end
end
text(5.45,0.05,'Probability of at least one measles case','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
ylim([0 0.1])


S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
M=1-p_zero_10;

Avg_Risk_County=zeros(length(S),1);
for ss=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
    if(sum(tf)>0)
        Avg_Risk_County(ss)=M(tf);
    end
end
Avg_Risk_County(Avg_Risk_County>0.1)=0.1;
NS=length(S);





CC_Cases=interp1(x_risk,C_Risk,Avg_Risk_County);
CC_Cases(isnan(Avg_Risk_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Risk_County)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.15,-0.15,1.2,1.2];

print(gcf,['Preliminary_County_Risk_Reduction_10.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% County
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=({'#ffffff';
     '#f7f4f9';
'#e7e1ef';
'#d4b9da';
'#c994c7';
'#df65b0';
'#e7298a';
'#ce1256';
'#980043';
'#67001f'; ...
'#000000'});


C_Risk=hex2rgb(C);
x_risk=linspace(0,0.1,11);


states = shaperead('usastatelo', 'UseGeoCoords', true);
figure('units','normalized','outerposition',[0 0.075 1 1]);

subplot("Position",[0.9 0.05 0.01 0.9])
avg_risk=linspace(0,0.1,1001);
CC_Risk=interp1(x_risk,C_Risk,avg_risk);

for jj=1:1000
    patch([0 0 1 1],[avg_risk(jj) avg_risk(jj+1) avg_risk(jj+1) avg_risk(jj)],CC_Risk(jj,:),'LineStyle','None'); hold on;
end

plot([0 1],[0 0],'k','LineWidth',2)
plot([0 1],0.1.*[1 1],'k','LineWidth',2)
plot([0 0],0.1.*[0 1],'k','LineWidth',2)
plot([1 1],0.1.*[0 1],'k','LineWidth',2)

text_v=[0:0.01:0.1];
for jj=1:length(text_v)
    if(jj==11)
        text(1.2,text_v(jj),['\geq ' num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    else
        text(1.2,text_v(jj),[num2str(text_v(jj))],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    end
end
text(5.45,0.05,'Probability of at least one measles case','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
ylim([0 0.1])


S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
M=1-p_zero;

Avg_Risk_County=zeros(length(S),1);
for ss=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(ss).GEOID);
    if(sum(tf)>0)
        Avg_Risk_County(ss)=M(tf);
    end
end
Avg_Risk_County(Avg_Risk_County>0.1)=0.1;
NS=length(S);





CC_Cases=interp1(x_risk,C_Risk,Avg_Risk_County);
CC_Cases(isnan(Avg_Risk_County),:)=repmat([0.5 0.5 0.5],sum(isnan(Avg_Risk_County)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Cases});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.15,-0.15,1.2,1.2];

print(gcf,['Preliminary_County_Risk_Baseline.png'],'-dpng','-r300');