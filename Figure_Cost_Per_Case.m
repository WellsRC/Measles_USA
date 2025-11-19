function Figure_Cost_Per_Case(Outbreak_Size)

NS=100;
load(['National_Reduction=' num2str(100*0) '_Year=' num2str(0) '.mat'],'County_Data_Vaccine_Reduction','Proportion_Size_Age_Unvaccinated','Proportion_Size_Age_Vaccinated');
[Total_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline]=Monte_Carlo_Fixed_Outbreak(Outbreak_Size,NS,County_Data_Vaccine_Reduction,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated);

[Total_Cost,Total_Outbreak_Response,Total_Medical,Total_Productivity_loss]=County_Level_Costs(County_Data_Vaccine_Reduction,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline);

Estimate_Cost_per_Case=median(Total_Cost./Total_Cases_County,2);


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
x_baseline=[75000 87700 100000];
C_Baseline=hex2rgb(Cb);

X_Label_Baseline=['Cost per case'];

prct_label=false;
monitary_label=true;

text_v=[75000 87700 100000];


inq_txt_baseline=zeros(size(text_v));
inq_txt_baseline(1)=-1;
inq_txt_baseline(end)=1;

Gen_Figure_Single_County(Cost_per_case_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S);

print(gcf,['Figure_3.png'],'-dpng','-r300');
end