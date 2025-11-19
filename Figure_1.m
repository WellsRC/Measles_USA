clear;
clc;

load('MMR_Immunity.mat','County_Data');
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
V_Baseline=NaN.*zeros(length(S),1);
U_State=unique(County_Data.State);
V_State=zeros(length(U_State),1);
for cc=1:length(S)
    tf=strcmp(County_Data.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        V_Baseline(cc)=(table2array(County_Data.Vaccine_Uptake(tf,1)).*table2array(County_Data.Population(tf,1))+(2/5).*table2array(County_Data.Vaccine_Uptake(tf,2)).*table2array(County_Data.Population(tf,2)))./(table2array(County_Data.Population(tf,1))+(2/5).*table2array(County_Data.Population(tf,2)));
    end
end

for ss=1:length(U_State)
    tf=strcmp(County_Data.State,U_State{ss});
    if(sum(tf)>0)

        V_temp=table2array(County_Data.Vaccine_Uptake(tf,1)).*table2array(County_Data.Population(tf,1)).*County_Data.Total_Population(tf)+(2/5).*table2array(County_Data.Vaccine_Uptake(tf,2)).*table2array(County_Data.Population(tf,2)).*County_Data.Total_Population(tf);
        P_temp=table2array(County_Data.Population(tf,1)).*County_Data.Total_Population(tf)+(2/5).*table2array(County_Data.Population(tf,2)).*County_Data.Total_Population(tf);

        V_State(ss)=sum(V_temp,1)./sum(P_temp,1);
    end   
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Plot Vaccine Uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 Cb=flip({'#ffffff';
     '#fff7ec';
'#fee8c8';
'#fdd49e';
'#fdbb84';
'#fc8d59';
'#ef6548';
'#d7301f';
'#b30000';
'#7f0000';});

 Measure_Baseline=V_Baseline;

x_baseline=linspace(0.65,0.95,length(Cb));
C_Baseline=hex2rgb(Cb);

X_Label_Baseline=['MMR coverage among children 0' char(8211) '6 years of age'];

prct_label=true;
monitary_label=false;

text_v=[0.65:0.1:0.95];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

Gen_Figure_Single_County(Measure_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S);

print(gcf,['Figure_1.png'],'-dpng','-r300');
