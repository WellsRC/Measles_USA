function Figure_Comparison_County(National_Reduction,Scenario,Age_0_to_6)
close all;
if(~Age_0_to_6)
    Age_Reduction=[true(1,1) false(1,4)];
else
    Age_Reduction=[true(1,2) false(1,3)];
end

[cases_baseline,hospital_baseline,cost_baseline,cost_per_case_baseline]=County_Outcome_Central_Measure(0,Scenario,Age_0_to_6);
[cases_reduction,hospital_reduction,cost_reduction,cost_per_case_reduction]=County_Outcome_Central_Measure(National_Reduction,Scenario,Age_0_to_6);
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

        Cost_Baseline(cc,:)=10.^3.*cost_baseline(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));
        Cost_Reduction(cc,:)=10.^3.*cost_reduction(tf)./(County_Data_Vaccine_Reduction.Total_Population(tf));

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Plot Vaccine Uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 Cb=flip({'#ffffff';
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

 Cd=flip({'#b2182b';
     '#ffffff';
     '#2166ac';});

Measure_Baseline=V_Baseline;
Measure_Reduction=V_Reduction;

x_baseline=linspace(0.5,1,11);
C_Baseline=hex2rgb(Cb);

x_diff=[0 National_Reduction, 0.02];
C_diff=hex2rgb(Cd);

X_Label_Baseline=['MMR uptake (0' char(8211) '6 years of age)'];
X_Label_Diff='Absolute reduction in MMR uptake';

prct_label=true;
monitary_label=false;

text_v=[0.5:0.1:1];
text_v_diff=[0 0.01 0.02];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;

inq_txt_diff=zeros(size(text_v_diff));

measure_increase=false;

Gen_Figure_Comparison(Measure_Baseline,Measure_Reduction,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,x_diff,C_diff,text_v_diff,inq_txt_diff,X_Label_Diff,S,measure_increase);

print(gcf,['Figure_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 Cb={'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d'};

 Cd=({'#ffffff';
'#fcc5c0';...
'#f768a1';
     '#ae017e';
'#49006a';});

Measure_Baseline=Case_Baseline;
Measure_Reduction=Case_Reduction;

x_baseline=[0.01 1:0.5:5];
C_Baseline=hex2rgb(Cb);

x_diff=[0:0.025:0.1];
C_diff=hex2rgb(Cd);

X_Label_Baseline='Measles cases per 10,000';
X_Label_Diff='Additional measles cases per 10,000';

prct_label=false;
monitary_label=false;

text_v=[0.01 1:5];
text_v_diff=[0:0.025:0.1];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

inq_txt_diff=zeros(size(text_v_diff));
inq_txt_diff(end)=1;

measure_increase=true;

Gen_Figure_Comparison(Measure_Baseline,Measure_Reduction,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,x_diff,C_diff,text_v_diff,inq_txt_diff,X_Label_Diff,S,measure_increase);

print(gcf,['Figure_Cases_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hospitalizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 Cb={'#ffffff';
     '#fff5f0';
'#fee0d2';
'#fcbba1';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#67000d'};

 Cd=({'#ffffff';
'#fcc5c0';
     '#dd3497';
'#49006a';});

Measure_Baseline=Hospital_Baseline;
Measure_Reduction=Hospital_Reduction;

x_baseline=[0.2 5:5:35];
C_Baseline=hex2rgb(Cb);

x_diff=[0.05 0.25:0.25:0.75];
C_diff=hex2rgb(Cd);

X_Label_Baseline='Hospitalizations per 1,000,000';
X_Label_Diff='Additional hospitalizations per 1,000,000';

prct_label=false;
monitary_label=false;

text_v=[0.2 5:5:35];
text_v_diff=[0.05 0.25:0.25:0.75];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

inq_txt_diff=zeros(size(text_v_diff));
inq_txt_diff(1)=-1;
inq_txt_diff(end)=1;

measure_increase=true;

Gen_Figure_Comparison(Measure_Baseline,Measure_Reduction,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,x_diff,C_diff,text_v_diff,inq_txt_diff,X_Label_Diff,S,measure_increase);


print(gcf,['Figure_Hospitalizations_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Costs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Cb={'#ffffff';
'#fee0d2';
'#fcbba1';
'#ef3b2c';
'#67000d'};

 Cd=({'#08306b'; ...
     '#ffffff';
'#fee0d2';
'#fcbba1';
'#fc9272';
'#fb6a4a';
'#ef3b2c';
'#cb181d';
'#a50f15';
'#67000d';});

Measure_Baseline=Cost_Baseline;
Measure_Reduction=Cost_Reduction;

x_baseline=[100 500 10^3 5.*10^3 10^4];
C_Baseline=hex2rgb(Cb);

x_diff=[-200 0 100:100:800];
C_diff=hex2rgb(Cd);

X_Label_Baseline='Cost of measles cases per 1,000';
X_Label_Diff='Additional cost of measles cases per 1,000';

prct_label=false;
monitary_label=true;

text_v=[10^2 5.*10^3 10^4];
text_v_diff=[-200 0 200:200:800];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

inq_txt_diff=zeros(size(text_v_diff));
inq_txt_diff(1)=-1;
inq_txt_diff(end)=1;

measure_increase=true;

Gen_Figure_Comparison(Measure_Baseline,Measure_Reduction,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,x_diff,C_diff,text_v_diff,inq_txt_diff,X_Label_Diff,S,measure_increase);


print(gcf,['Figure_Cost_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost_per_case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 Cb=flip({'#543005';
'#8c510a';
'#bf812d';
'#dfc27d';
'#f6e8c3';
'#f5f5f5';
'#c7eae5';
'#80cdc1';
'#35978f';
'#01665e';
'#003c30'});

 Cd=flip({'#b2182b';
'#d6604d';
'#f4a582';
'#fddbc7';
'#ffffff';
'#d1e5f0';
'#92c5de';
'#4393c3';
'#2166ac';});

Measure_Baseline=Cost_per_case_Baseline;
Measure_Reduction=Cost_per_case_Reduction;

x_baseline=[42500 45000 47500 50000 52500 55000 57500 60000 62500 65000 67500];
C_Baseline=hex2rgb(Cb);

x_diff=[-500:125:500];
C_diff=hex2rgb(Cd);

X_Label_Baseline='Cost per measles case';
X_Label_Diff='Additional cost per measles case';

prct_label=false;
monitary_label=true;

text_v=[42500 55000 67500];
text_v_diff=[-500:250:500];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

inq_txt_diff=zeros(size(text_v_diff));
inq_txt_diff(1)=-1;
inq_txt_diff(end)=1;

measure_increase=true;

Gen_Figure_Comparison(Measure_Baseline,Measure_Reduction,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,x_diff,C_diff,text_v_diff,inq_txt_diff,X_Label_Diff,S,measure_increase);

print(gcf,['Figure_Cost_per_case_Vaccine_Uptake_Reduction=' num2str(National_Reduction*100) '.png'],'-dpng','-r300');
end