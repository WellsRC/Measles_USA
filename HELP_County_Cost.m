function HELP_County_Cost()
close all;
Scenario='Sample_2025';
Age_0_to_6=true;

[cases_baseline,hospital_baseline,cost_baseline,cost_per_case_baseline]=County_Outcome_Central_Measure(0,Scenario,Age_0_to_6);
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);

load(['National_Reduction=' num2str(100*0) '_Ages_0_to_6.mat'],'County_Data_Vaccine_Reduction')

V_Baseline=NaN.*zeros(length(S),1);

Case_Baseline=NaN.*zeros(length(S),1);

Hospital_Baseline=NaN.*zeros(length(S),1);

Cost_Baseline=NaN.*zeros(length(S),1);

Cost_per_case_Baseline=NaN.*zeros(length(S),1);

for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        Cost_per_case_Baseline(cc)=cost_per_case_baseline(tf);
        V_Baseline(cc)=(table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,1)).*table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(tf,2)).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)))./(table2array(County_Data_Vaccine_Reduction.Population(tf,1))+(2/5).*table2array(County_Data_Vaccine_Reduction.Population(tf,2)));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Plot Vaccine Uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 Cb={'#003c30'; ...
     '#ffffff';
'#543005'};

 Measure_Baseline=Cost_per_case_Baseline;

x_baseline=[79000 91500 104000];
C_Baseline=hex2rgb(Cb);

X_Label_Baseline=['Average cost per case'];

prct_label=false;
monitary_label=true;

text_v=[79000 91500 104000];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

Gen_Figure_Single_County(Measure_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S);

print(gcf,['HELP_County_Cost_per_Case.png'],'-dpng','-r300');

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

X_Label_Baseline=['MMR uptake among children 0' char(8211) '6 years of age'];

prct_label=true;
monitary_label=false;

text_v=[0.65:0.1:0.95];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

Gen_Figure_Single_County(Measure_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S);

print(gcf,['HELP_MMR_Uptake_0_to_6.png'],'-dpng','-r300');
end