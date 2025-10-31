clear;
[Importation_Cases_County] = Case_Importation_Sample('Sample_2025',2.5.*10^3);
pd=fitdist(sum(Importation_Cases_County,1)','Kernel','Support','positive');
xm=fmincon(@(x)-pdf(pd,x),mean(sum(Importation_Cases_County,1)),[],[],[],[],0,max(sum(Importation_Cases_County,1)));

Year=[2023:2025];