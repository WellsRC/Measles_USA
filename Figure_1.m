clear;
clc;

Outbreak_Size=50;
NS=2500;
load(['National_Reduction=' num2str(0) '_Year=' num2str(0) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated','Proportion_Age_Unvaccinated_Uninsured','Proportion_Age_Unvaccinated_Public','Proportion_Age_Unvaccinated_Private','Proportion_Age_Vaccinated_Uninsured','Proportion_Age_Vaccinated_Public','Proportion_Age_Vaccinated_Private');
[Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts]=Monte_Carlo_Fixed_Outbreak(Outbreak_Size,NS,County_Data_Vaccine_Reduction,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,Proportion_Age_Unvaccinated_Uninsured,Proportion_Age_Unvaccinated_Public,Proportion_Age_Unvaccinated_Private,Proportion_Age_Vaccinated_Uninsured,Proportion_Age_Vaccinated_Public,Proportion_Age_Vaccinated_Private);

[Total_Cost,~,~,~]=County_Level_Costs(County_Data_Vaccine_Reduction,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County,Total_Contacts,Unvaccinated_Contacts);

Estimate_Cost_per_Case=median(Total_Cost./Total_Cases_County,2);

[r,p]=corr(County_Data_Vaccine_Reduction.Total_Immunity(:),Estimate_Cost_per_Case(:),'Type','Spearman')

std(Estimate_Cost_per_Case)
prctile(Estimate_Cost_per_Case,[50 25 75])

 S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);

Cost_per_case_Baseline=NaN.*zeros(length(S),1);
for cc=1:length(S)
    tf=strcmp(County_Data_Vaccine_Reduction.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        Cost_per_case_Baseline(cc)=Estimate_Cost_per_Case(tf);
    end
end


Cb={'#003c30'; ...
     '#ffffff';
'#543005'};
x_baseline=[75000 87000 101000];
C_Baseline=hex2rgb(Cb);

X_Label_Baseline=['Cost per case'];

prct_label=false;
monitary_label=true;

text_v={'$75,000','$87,000','$101,000'};


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

f=Gen_Figure_Single_County(Cost_per_case_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S,0);

theme(f, "light"); 
print(gcf,['Figure_1.png'],'-dpng','-r300');
print(gcf,['Figure_1.tif'],'-dtiff','-r300');
