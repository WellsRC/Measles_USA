function Gen_Figure_Single_County(Measure_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S,label_plot)


states = shaperead('usastatelo', 'UseGeoCoords', true);

figure('units','normalized','outerposition',[0.15 0.075 0.3 0.43]);

subplot("Position",[0.085 0.135 0.83 0.025])
dv=linspace(x_baseline(1),x_baseline(end),1001);
CC_Baseline=interp1(x_baseline,C_Baseline,dv);

for jj=1:1000
    patch([dv(jj) dv(jj+1) dv(jj+1) dv(jj)],[0 0 1 1],CC_Baseline(jj,:),'LineStyle','None'); hold on;
end

plot([x_baseline(1) x_baseline(1)],[0 1],'k','LineWidth',1)
plot([x_baseline(end) x_baseline(end)],[0 1],'k','LineWidth',1)
plot([x_baseline(1) x_baseline(end)],[0 0],'k','LineWidth',1)
plot([x_baseline(1) x_baseline(end)],[1 1],'k','LineWidth',1)

if(prct_label)
    for jj=1:length(text_v)
        if(inq_txt_baseline(jj)<0)
            text(text_v(jj),-1.2,['\leq ' num2str(100.*abs(text_v(jj))) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        elseif(inq_txt_baseline(jj)>0)
            text(text_v(jj),-1.2,['\geq ' num2str(100.*abs(text_v(jj))) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        else
            text(text_v(jj),-1.2,[num2str(100.*abs(text_v(jj))) '%'],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        end
    end
elseif(monitary_label)
    for jj=1:length(text_v)
        if(text_v(jj)<0)
            sgn='-';
        else
            sgn=[];
        end
        if(inq_txt_baseline(jj)<0)
            text(text_v(jj),-1.2,['\leq' sgn '$' num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        elseif(inq_txt_baseline(jj)>0)
            text(text_v(jj),-1.2,['\geq' sgn '$' num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        else
            text(text_v(jj),-1.2,[sgn '$' num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        end
    end
else
    for jj=1:length(text_v)
        if(inq_txt_baseline(jj)<0)
            text(text_v(jj),-1.2,['\leq ' num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        elseif(inq_txt_baseline(jj)>0)
            text(text_v(jj),-1.2,['\geq ' num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        else
            text(text_v(jj),-1.2,[num2str(abs(text_v(jj)))],'Fontsize',14,'VerticalAlignment','middle','HorizontalAlignment','center');
        end
    end
end
text(0.5,-3.9,X_Label_Baseline,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle','Units','normalized');
if(label_plot>0)
    text(-0.099,32.28,char(64+label_plot),'Fontsize',32,'Units','normalized')
end
axis off;
xlim([x_baseline(1) x_baseline(end)])

Measure_Baseline(Measure_Baseline<x_baseline(1))=x_baseline(1);
Measure_Baseline(Measure_Baseline>x_baseline(end))=x_baseline(end);

NS=length(S);

CC_Baseline=interp1(x_baseline,C_Baseline,Measure_Baseline(:));
CC_Baseline(isnan(Measure_Baseline(:)),:)=repmat([0.5 0.5 0.5],sum(isnan(Measure_Baseline(:))),1);
CM=makesymbolspec('Polygon',{'INDEX',[1 NS],'FaceColor',CC_Baseline});

ax1=usamap('conus');
 framem off; gridm off; mlabel off; plabel off;
ax1.Position=[-0.3,-2,0.6,0.6];
geoshow(ax1,S,'SymbolSpec',CM,'LineStyle','None'); hold on;
geoshow(ax1, states,'Facecolor','none','LineWidth',1.5); 

ax1.Position=[-0.135, -0.11 ,1.3,1.3];



end