function Figure_Vaccine_Uptake_Reductions(reduction_uptake)
clear;
close all;
Vaccine='MMR';
L_F=zeros(2^7,4);
L_S=zeros(2^7,4);
for Spatial_Valiation=1:4
    load(['Refine_' Vaccine '_Set_Spatial_Validation=' num2str(Spatial_Valiation) '.mat']);
    L_F(:,Spatial_Valiation)=L_fit;
    L_S(:,Spatial_Valiation)=L_spatial_val;
end

L_T=sum(L_F,2)+sum(L_S,2);
    
[~,Model_Num]=max(L_T);

clearvars -except Model_Num Vaccine Age_Class


[Age_0_to_4_county]=Age_Adjustment_Factor(Vaccine,Model_Num,'Age_0_to_4');
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
GEOID_0_to_4=Age_0_to_4_county.GEOID;
Vaccine_Uptake_0_to_4=NaN.*zeros(length(S),1);

for ss=1:length(S)
    tf=strcmp(GEOID_0_to_4,S(ss).GEOID);
    if(sum(tf)>0)
        Vaccine_Uptake_0_to_4(ss)=Age_0_to_4_county.Vaccine_Uptake(tf);
    end
end

[Age_5_to_9_county]=Age_Adjustment_Factor(Vaccine,Model_Num,'Age_5_to_9');
GEOID_5_to_9=Age_5_to_9_county.GEOID;
Vaccine_Uptake_5_to_9=NaN.*zeros(length(S),1);

for ss=1:length(S)
    tf=strcmp(GEOID_5_to_9,S(ss).GEOID);
    if(sum(tf)>0)
        Vaccine_Uptake_5_to_9(ss)=Age_5_to_9_county.Vaccine_Uptake(tf);
    end
end

[Age_10_to_14_county]=Age_Adjustment_Factor(Vaccine,Model_Num,'Age_10_to_14');
GEOID_10_to_14=Age_10_to_14_county.GEOID;
Vaccine_Uptake_10_to_14=NaN.*zeros(length(S),1);

for ss=1:length(S)
    tf=strcmp(GEOID_10_to_14,S(ss).GEOID);
    if(sum(tf)>0)
        Vaccine_Uptake_10_to_14(ss)=Age_10_to_14_county.Vaccine_Uptake(tf);
    end
end

[Age_15_to_19_county]=Age_Adjustment_Factor(Vaccine,Model_Num,'Age_15_to_19');
GEOID_15_to_19=Age_15_to_19_county.GEOID;
Vaccine_Uptake_15_to_19=NaN.*zeros(length(S),1);

for ss=1:length(S)
    tf=strcmp(GEOID_15_to_19,S(ss).GEOID);
    if(sum(tf)>0)
        Vaccine_Uptake_15_to_19(ss)=Age_15_to_19_county.Vaccine_Uptake(tf);
    end
end

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


states = shaperead('usastatelo', 'UseGeoCoords', true);
figure('units','normalized','outerposition',[0 0.075 1 1]);

subplot("Position",[0.9 0.05 0.01 0.9])
dv=linspace(0.05,1,1001);
CC_Vac=interp1(x_vac,C_vac,dv);

for jj=1:1000
    patch([0 0 1 1],100.*[dv(jj) dv(jj+1) dv(jj+1) dv(jj)],CC_Vac(jj,:),'LineStyle','None'); hold on;
end

plot([0 1],[5 5],'k','LineWidth',2)
plot([0 1],[100 100],'k','LineWidth',2)
plot([0 0],[5 100],'k','LineWidth',2)
plot([1 1],[5 100],'k','LineWidth',2)

text_v=[5:5:100];
for jj=1:length(text_v)
    if(jj==1)
        text(1.2,text_v(jj),['\leq ' num2str(text_v(jj)) '%'],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    else
        text(1.2,text_v(jj),[num2str(text_v(jj)) '%'],'Fontsize',20,'VerticalAlignment','middle','HorizontalAlignment','left');
    end
end
text(5.45,52.5,'MMR vaccine uptake','Rotation',270,'FontSize',32,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
ylim([5 100])


Vaccine_Uptake_0_to_4(Vaccine_Uptake_0_to_4<0.05)=0.05;
Vaccine_Uptake_5_to_9(Vaccine_Uptake_5_to_9<0.05)=0.05;
Vaccine_Uptake_10_to_14(Vaccine_Uptake_10_to_14<0.05)=0.05;
Vaccine_Uptake_15_to_19(Vaccine_Uptake_15_to_19<0.05)=0.05;


NS=length(S);





CC_Vac=interp1(x_vac,C_vac,Vaccine_Uptake_0_to_4);
CC_Vac(isnan(Vaccine_Uptake_0_to_4),:)=repmat([0.5 0.5 0.5],sum(isnan(Vaccine_Uptake_0_to_4)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,Vaccine_Uptake_5_to_9);
CC_Vac(isnan(Vaccine_Uptake_5_to_9),:)=repmat([0.5 0.5 0.5],sum(isnan(Vaccine_Uptake_5_to_9)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax2=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax2.Position=[1.7,0.4,0.6,0.6];
geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,Vaccine_Uptake_10_to_14);
CC_Vac(isnan(Vaccine_Uptake_10_to_14),:)=repmat([0.5 0.5 0.5],sum(isnan(Vaccine_Uptake_10_to_14)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax3=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax3.Position=[-0.3,-0.1,0.6,0.6];

geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,Vaccine_Uptake_15_to_19);
CC_Vac(isnan(Vaccine_Uptake_15_to_19),:)=repmat([0.5 0.5 0.5],sum(isnan(Vaccine_Uptake_15_to_19)),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax4=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax4.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax4,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax4, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.075,0.425,0.6,0.6];
ax2.Position=[0.35,0.425,0.6,0.6];
ax3.Position=[-0.075,-0.075,0.6,0.6];
ax4.Position=[0.35,-0.075,0.6,0.6];

text(-0.65,1.75,'A','Units','normalized','Fontsize',40)
text(-0.42,1.75,'Ages 0 to 4','Units','normalized','Fontsize',32)

text(0.135,1.75,'B','Units','normalized','Fontsize',40)
text(0.35,1.75,'Ages 5 to 9','Units','normalized','Fontsize',32)

text(-0.65,0.915,'C','Units','normalized','Fontsize',40)
text(-0.42,0.915,'Ages 10 to 14','Units','normalized','Fontsize',32)

text(0.135,0.915,'D','Units','normalized','Fontsize',40)
text(0.35,0.915,'Ages 15 to 19','Units','normalized','Fontsize',32)

end