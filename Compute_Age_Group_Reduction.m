function Age_Group_Reduction = Compute_Age_Group_Reduction(National_Annual_Reduction,Year_Reduced)
load(['National_Reduction=0_Year=0.mat'],'County_Data_Vaccine_Reduction')

Pop_Single=readtable('Child population by single age.xlsx');
Pop_Single=Pop_Single(strcmp(Pop_Single.LocationType,'State') & Pop_Single.TimeFrame==2023 & Pop_Single.SingleAge<15,:);
P_Age=zeros(1,15);
for ii=1:3
    for jj=1:5
        P_Age(jj+5.*(ii-1))=sum(Pop_Single.Data(Pop_Single.SingleAge==5.*(ii-1)+(jj-1)));
    end
end

V=County_Data_Vaccine_Reduction.Total_Population'*(table2array(County_Data_Vaccine_Reduction.Population(:,1:5)).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1:5)));
VC=V./(County_Data_Vaccine_Reduction.Total_Population'*(table2array(County_Data_Vaccine_Reduction.Population(:,1:5))));

V_NAT_0_4=[0 0.654 0.906 0.931 0.927]; % Age 1: 13 month 2021 bith-year Age 2: 24 month 2021 bith-year Age 3: 35 month 2021 bith-year Age 4: 35 month 2020 bith-year for national value

V0_Age=[V_NAT_0_4.*P_Age(1:5) VC(2).*P_Age(6:10) VC(3).*P_Age(11:15)];

VC0_Age=V0_Age./P_Age;


P_0_to_6=P_Age(1:7);
VB1=[V0_Age(1:7)];


VBT=sum(VB1);

New_Vac=VC0_Age;

max_A=7+(Year_Reduced-1);
r=lsqnonlin(@(x) sum((VC0_Age(1:7).*[0 1./(1+exp(-x./(max_A-[1:6])))]).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-National_Annual_Reduction.*Year_Reduced),14,0,100);

New_Vac(1:(max_A+1))=VC0_Age(1:(max_A+1)).*[0 1./(1+exp(-r./(max_A-[1:max_A])))];

Age_Group_Reduction=zeros(1,3);
for ii=1:3
    Age_Group_Reduction(ii)=sum(VC0_Age(5.*(ii-1)+[1:5]).*P_Age(5.*(ii-1)+[1:5]))./sum(P_Age(5.*(ii-1)+[1:5]))-sum(New_Vac(5.*(ii-1)+[1:5]).*P_Age(5.*(ii-1)+[1:5]))./sum(P_Age(5.*(ii-1)+[1:5]));
end

end