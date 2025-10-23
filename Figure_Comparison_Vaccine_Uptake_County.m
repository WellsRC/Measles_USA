function Figure_Comparison_Vaccine_Uptake_County(National_Reduction,Age_Reduction)
close all;
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);

FN_Age_Class={'Age_0_to_4','_Age_5_to_9','_Age_10_to_14','_Age_15_to_19','_Age_20_to_24'};
FN_Age_Class=FN_Age_Class(Age_Reduction);
load(['National_Reduction=' num2str(100*0) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')

V_Baseline=NaN.*zeros(length(S),sum(Age_Reduction));
V_Reduction=NaN.*zeros(length(S),sum(Age_Reduction));

for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        V_Baseline(cc,:)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
    end
end

load(['National_Reduction=' num2str(100*National_Reduction) '_' FN_Age_Class{:} '.mat'],'County_Data_Vaccine_Reduction')
for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        V_Reduction(cc,:)=table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,Age_Reduction));
    end
end

V_Diff=V_Baseline-V_Reduction;


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



 C={'#ffffff';
     '#fff7f3';
'#fde0dd';
'#fcc5c0';
'#fa9fb5';
'#f768a1';
'#dd3497';
'#ae017e';
'#7a0177';
'#49006a'};


C_diff=hex2rgb(C);
x_diff=linspace(0,ceil(100.*(max(V_Diff(:))))./100,size(C_diff,1));

states = shaperead('usastatelo', 'UseGeoCoords', true);

figure('units','normalized','outerposition',[0.15 0.075 0.6 0.9]);

subplot("Position",[0.0225 0.065 0.28 0.015])
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
text(52.5,-2.9,'MMR vaccine uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
xlim([5 100])

subplot("Position",[0.36 0.065 0.28 0.015])
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
text(52.5,-2.9,'MMR vaccine uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
xlim([5 100])


subplot("Position",[0.69 0.065 0.28 0.015])
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
text(50.*max(x_diff),-2.9,'Difference in uptake','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle')
axis off;
xlim([0 100.*max(x_diff)])

V_Baseline(V_Baseline<0.05)=0.05;
V_Reduction(V_Reduction<0.05)=0.05;




NS=length(S);




Age_Count=1;
CC_Vac=interp1(x_vac,C_vac,V_Baseline(:,Age_Count));
CC_Vac(isnan(V_Baseline(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Baseline(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,V_Reduction(:,Age_Count));
CC_Vac(isnan(V_Reduction(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Reduction(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});
ax2=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax2.Position=[1.7,0.4,0.6,0.6];
geoshow(ax2,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax2, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,V_Diff(:,Age_Count));
CC_Diff(isnan(V_Diff(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Diff(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax3=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax3.Position=[-0.3,-0.1,0.6,0.6];

geoshow(ax3,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax3, states,'Facecolor','none','LineWidth',1.5); 


Age_Count=2;
CC_Vac=interp1(x_vac,C_vac,V_Baseline(:,Age_Count));
CC_Vac(isnan(V_Baseline(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Baseline(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax4=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax4.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax4,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax4, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,V_Reduction(:,Age_Count));
CC_Vac(isnan(V_Reduction(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Reduction(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax5=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax5.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax5,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax5, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,V_Diff(:,Age_Count));
CC_Diff(isnan(V_Diff(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Diff(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax6=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax6.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax6,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax6, states,'Facecolor','none','LineWidth',1.5); 


Age_Count=3;
CC_Vac=interp1(x_vac,C_vac,V_Baseline(:,Age_Count));
CC_Vac(isnan(V_Baseline(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Baseline(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax7=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax7.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax7,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax7, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,V_Reduction(:,Age_Count));
CC_Vac(isnan(V_Reduction(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Reduction(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax8=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax8.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax8,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax8, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,V_Diff(:,Age_Count));
CC_Diff(isnan(V_Diff(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Diff(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax9=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax9.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax9,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax9, states,'Facecolor','none','LineWidth',1.5); 

Age_Count=4;
CC_Vac=interp1(x_vac,C_vac,V_Baseline(:,Age_Count));
CC_Vac(isnan(V_Baseline(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Baseline(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax10=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax10.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax10,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax10, states,'Facecolor','none','LineWidth',1.5); 


CC_Vac=interp1(x_vac,C_vac,V_Reduction(:,Age_Count));
CC_Vac(isnan(V_Reduction(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Reduction(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Vac});

ax11=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax11.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax11,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax11, states,'Facecolor','none','LineWidth',1.5); 


CC_Diff=interp1(x_diff,C_diff,V_Diff(:,Age_Count));
CC_Diff(isnan(V_Diff(:,Age_Count)),:)=repmat([0.5 0.5 0.5],sum(isnan(V_Diff(:,Age_Count))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Diff});

ax12=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;

ax12.Position=[1.7,-0.1,0.6,0.6];
geoshow(ax12,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax12, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.01, 0.66,0.4,0.4];
ax2.Position=[0.32,0.66,0.4,0.4];
ax3.Position=[0.63, 0.66,0.4,0.4];

ax4.Position=[-0.01, 0.44,0.4,0.4];
ax5.Position=[0.32,0.44,0.4,0.4];
ax6.Position=[0.63, 0.44,0.4,0.4];

ax7.Position=[-0.01, 0.2,0.4,0.4];
ax8.Position=[0.32,0.22,0.4,0.4];
ax9.Position=[0.63, 0.22,0.4,0.4];

ax10.Position=[-0.01, 0.0,0.4,0.4];
ax11.Position=[0.32,0.0,0.4,0.4];
ax12.Position=[0.63, 0.0,0.4,0.4];

end

