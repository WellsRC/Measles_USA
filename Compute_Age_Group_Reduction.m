function Age_Group_Reduction = Compute_Age_Group_Reduction(National_Annual_Reduction,Year_Reduced)
load(['National_Reduction=0_Year=0.mat'],'County_Data_Vaccine_Reduction')

Pop_Single=readtable('Child population by single age.xlsx');
Pop_Single=Pop_Single(strcmp(Pop_Single.LocationType,'State') & Pop_Single.TimeFrame==2023 & Pop_Single.SingleAge<15,:);
P_Age=zeros(3,5);
for ii=1:3
    for jj=1:5
        P_Age(ii,jj)=sum(Pop_Single.Data(Pop_Single.SingleAge==5.*(ii-1)+(jj-1)));
    end
end

V=County_Data_Vaccine_Reduction.Total_Population'*(table2array(County_Data_Vaccine_Reduction.Population(:,1:5)).*table2array(County_Data_Vaccine_Reduction.Vaccine_Uptake(:,1:5)));
VC=V./(County_Data_Vaccine_Reduction.Total_Population'*(table2array(County_Data_Vaccine_Reduction.Population(:,1:5))));

V_NAT_0_to_4=[0 0.654 0.906 0.931 0.927]; % Age 1: 13 month 2021 bith-year Age 2: 24 month 2021 bith-year Age 3: 35 month 2021 bith-year Age 4: 35 month 2020 bith-year for national value

V0_Age=[V_NAT_0_to_4.*P_Age(1,:);VC(2).*P_Age(2,:);VC(3).*P_Age(3,:)];

VC0_Age=V0_Age./P_Age;


P_0_to_6=[P_Age(1,1:5) P_Age(2,1:2)];
VB1=[V0_Age(1,:) V0_Age(2,1:2)];


VBT=sum(VB1);

VC0=VB1./P_0_to_6;

lb=[0 VC0(1:6)];
ub=[VC0(1:7)];

r=lsqnonlin(@(x) sum((lb+(ub-lb).*x).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-National_Annual_Reduction),0.9,0,1);

V1=[(lb+(ub-lb).*r).*P_0_to_6];

VC1=[V1./P_0_to_6 VC0(7:end)];


VC1_Age=VC0_Age;
VC1_Age(1,:)=[VC1(1:5)];
VC1_Age(2,1:3)=VC1(6:end);


lb=[0 VC1(1:6)];
ub=[VC1(1:7)];


r=lsqnonlin(@(x) sum((lb+(ub-lb).*x).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-2.*National_Annual_Reduction),0.9,0,1);

V2=[(lb+(ub-lb).*r).*P_0_to_6];
VC2=[V2./P_0_to_6 VC1(7:end)];

VC2_Age=VC0_Age;
VC2_Age(1,:)=[VC2(1:5)];
VC2_Age(2,1:4)=VC2(6:end);

lb=[0 VC2(1:6)];
ub=[VC2(1:7)];
r=lsqnonlin(@(x) sum((lb+(ub-lb).*x).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-3.*National_Annual_Reduction),0.9,0,1);


V3=[(lb+(ub-lb).*r).*P_0_to_6];


VC3=[V3./P_0_to_6 VC2(7:end)];

VC3_Age=VC0_Age;
VC3_Age(1,:)=[VC3(1:5)];
VC3_Age(2,1:5)=VC3(6:end);

lb=[0 VC3(1:6)];
ub=[VC3(1:7)];
r=lsqnonlin(@(x) sum((lb+(ub-lb).*x).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-4.*National_Annual_Reduction),0.9,0,1);


V4=[(lb+(ub-lb).*r).*P_0_to_6];

VC4=[V4./P_0_to_6 VC3(7:end)];

VC4_Age=VC0_Age;
VC4_Age(1,:)=[VC4(1:5)];
VC4_Age(2,:)=VC4(6:10);
VC4_Age(3,1)=VC4(11:end);

lb=[0 VC4(1:6)];
ub=[VC4(1:7)];
r=lsqnonlin(@(x) sum((lb+(ub-lb).*x).*P_0_to_6)./sum(P_0_to_6)-(VBT./sum(P_0_to_6)-5.*National_Annual_Reduction),0.9,0,1);


V5=[(lb+(ub-lb).*r).*P_0_to_6];

VC5=[V5./P_0_to_6 VC4(7:end)];

VC5_Age=VC0_Age;
VC5_Age(1,:)=[VC5(1:5)];
VC5_Age(2,:)=VC5(6:10);
VC5_Age(3,1:2)=VC5(11:end);

Reductions_Age_Groups=zeros(5,3);

for aa=1:3
    Pt=[P_Age(aa,:)];
test0=[VC0_Age(aa,:)];
test1=[VC1_Age(aa,:)];
test2=[VC2_Age(aa,:)];
test3=[VC3_Age(aa,:)];
test4=[VC4_Age(aa,:)];
test5=[VC5_Age(aa,:)];

Reductions_Age_Groups(1,aa)=(Pt*test0(:)./sum(Pt))-(Pt*test1(:)./sum(Pt));
Reductions_Age_Groups(2,aa)=(Pt*test0(:)./sum(Pt))-(Pt*test2(:)./sum(Pt));
Reductions_Age_Groups(3,aa)=(Pt*test0(:)./sum(Pt))-(Pt*test3(:)./sum(Pt));
Reductions_Age_Groups(4,aa)=(Pt*test0(:)./sum(Pt))-(Pt*test4(:)./sum(Pt));
Reductions_Age_Groups(5,aa)=(Pt*test0(:)./sum(Pt))-(Pt*test5(:)./sum(Pt));

end


Age_Group_Reduction=Reductions_Age_Groups(Year_Reduced,:);
end